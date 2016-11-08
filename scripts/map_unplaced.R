setwd("~/Dropbox/pmdvlab/do_fix/")

library(ggplot2)
library(plyr)
library(qtl2scan)
library(qtl2geno)
library(GenomicRanges)
devtools::load_all("~/Dropbox/pmdvlab/mouser/")

source("~/lib/util/R/wgs.R")
source("~/lib/util/R/qtl2.R")
source("~/Dropbox/pmdvlab/R_multiparental/R/ggqtlplot.R")

covg <- read.coverage("wgs/unplaced/all.unplaced.w10kb.bed", has.MQ0 = TRUE)
covg$iid <- gsub("_", ".", covg$iid)

covg2 <- ddply(covg, .(chr, iid), summarise,
			   lo = min(depth, na.rm = TRUE), hi = max(depth, na.rm = TRUE),
			   med = median(depth, na.rm = TRUE), mu = mean(depth, na.rm = TRUE),
			   lq = quantile(depth, 0.25, na.rm = TRUE), uq = quantile(depth, 0.74, na.rm = TRUE))

covg3 <- ddply(covg, .(chr, start, end), summarise,
			   med = median(depth, na.rm = TRUE), lq = quantile(depth, 0.25, na.rm = TRUE),
			   uq = quantile(depth, 0.75, na.rm = TRUE))
covg3 <- ddply(covg3, .(chr), transform, idx = seq_along(start))

pdf("wgs/unplaced/depth.bysample.pdf", width = 8, height = 6)
d_ply(covg2, .(chr), function(df) {
	p <- ggplot(df) +
		geom_pointrange(aes(x = reorder(iid, med), y = med, ymin = lq, ymax = uq)) +
		geom_hline(yintercept = 2.0, lty = "dashed", colour = "blue") +
		scale_y_continuous("normalized read depth\n", trans = "log2") +
		scale_x_discrete("\nsamples (ordered by coverage)") +
		facet_grid(. ~ chr) +
		theme_bw() +
		theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
	print(p)
}, .progress = "text")
dev.off()

pdf("wgs/unplaced/depth.bycontig.pdf", width = 8, height = 6)
d_ply(covg3, .(chr), function(df) {
	p <- ggplot(df) +
		geom_ribbon(aes(x = idx, ymin = lq, ymax = uq), fill = "grey") +
		geom_line(aes(x = idx, y = med)) +
		scale_y_continuous("normalized read depth\n", trans = "log2") +
		facet_wrap(~ chr) +
		theme_bw()
	print(p)
}, .progress = "text")
dev.off()

## load samples
samples <- openxlsx::read.xlsx("wgs/wgs_samples.xlsx", 1)
rownames(samples) <- samples$accession

## load genotype probs, qtl2 format
geno <- readRDS("wgs/genoprobs.qtl2.rds")
K <- readRDS("wgs/kinship.rds")
all(dimnames(K)[[1]] == dimnames(geno)[[1]])
all(rownames(K) %in% samples$iid)

## load SNPs info (mm10 coords)
snps <- readRDS("~/db/arrays/megamuga/snps.megamuga.mm10.rds")
snps <- snps[ allmarkers(geno), ]
snps2 <- data.frame(marker = snps$marker, chr = gsub("^chr","", snps$chr),
					pos = snps$pos.mm10/1e6, cM = snps$cM,
					row.names = as.character(snps$marker))

## for mapping, phenotype = median copy number across a segment
## (per plots above, this is not too bad a proxy for the real copy number)
pheno <- reshape2::acast(covg2, iid ~ chr, value.var = "med")
rownames(pheno) <- as.character(pheno$iid)
pheno <- matrix(pheno[ ,c("med") ], ncol = 1, dimnames = list(rownames(pheno)))

## run genome scans and add SNP info
rez <- scan1(geno, pheno)
rez2 <- tidy.scan1(rez)
rez2 <- merge(rez2, snps[ ,c("marker","pos") ])

## extract peaks; estimate strain means
peaks <- qtl2scan::find_peaks(rez, threshold = 10, drop = 5)
betas <- ddply(peaks, .(lodindex, chr), effects.at.peaks, geno = geno, pheno = pheno)
betas$trait <- colnames(rez$lod)[ betas$lodindex ]
betas$label <- locus.label(betas$chr, betas$cM, "cM")

## make a bunch of plots

pdf("wgs/unplaced/mapping.lodcurves.pdf", width = 10, height = 5)
d_ply(rez2, .(trait), function(df) {
	toplot <- merge(df, snps[ ,c("marker","pos") ])
	p <- ggqtlplot(toplot, marker.rug = FALSE) +
		ggtitle(as.character(df$trait[1]))
	print(p)
}, .progress = "text")
dev.off()

pdf("wgs/unplaced/mapping.straineff.pdf", width = 8, height = 10)
d_ply(betas, .(trait), function(df) {
	
	the.trait <- as.character(df$trait)[1]
	toplot <- merge(subset(rez2, trait == the.trait), snps[ ,c("marker","pos") ])
	p1 <- ggqtlplot(toplot, marker.rug = FALSE) +
		ggtitle(the.trait)
	
	p2 <- ggplot(df) +
		geom_pointrange(aes(x = strain, y = beta, ymin = beta-2*se, ymax = beta+2*se, colour = strain)) +
		geom_hline(yintercept = 2.0, lty = "dashed", colour = "black") +
		scale_y_continuous("copy number\n") +
		scale_color_CC() +
		facet_grid(. ~ label)
	
	print( cowplot::plot_grid(p1, p2, nrow = 2) )
	
}, .progress = "text")
dev.off()