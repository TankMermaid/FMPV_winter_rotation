# qtl2.R
# utilities for doing genome scans with from R/qtl2 and plotting results

#' Extract marker names on all chromosomes.
allmarkers <- function(gprobs, ...) {
	do.call("c", lapply(gprobs$map, names))
}

#' Extract genotype probabilities at a single marker.
grab.marker <- function(gprobs, mk, ...) {
	probs <- gprobs$probs[[1]][ ,,mk, drop = FALSE ]
	map <- gprobs$map[[1]][mk]
	gprobs$probs[[1]] <- probs
	gprobs$map[[1]] <- map
	return(gprobs)
}

#' Estimate allele effects at given QTL peaks from qtl2scan::find_peaks()
#' @details This function is intended for use with plyr: I expect each input to correspond to a single phenotype,
#' 	for example: ddply(peaks, .(lodindex), effects.at.peaks, ...)
effects.at.peaks <- function(peaks, geno, pheno, se = TRUE, ...) {
	
	chrom <- peaks$chr[1]
	this.geno <- subset(geno, chr = chrom)
	map <- sort(this.geno$map[[1]])
	ii <- findInterval(peaks$pos, map)
	mk <- names(map)[ii]
	
	coefs <- scan1blup(grab.marker(this.geno, mk), pheno[ ,peaks$lodindex[1] ], se = se)
	tidy.scan1coef(coefs)
	
}

#' Tidy up coefficient estimates from qtl2scan::scan1blup()
tidy.scan1coef <- function(x, ...) {
	
	coefs <- reshape2::melt(x$coef)
	colnames(coefs) <- c("marker","strain","beta")
	if (!is.null(x$SE)) {
		se <- reshape2::melt(x$SE)
		colnames(se) <- c("marker","strain","se")
		coefs <- cbind(coefs, se = se[,"se"])
	}
	
	map <- data.frame(marker = names(x$map), cM = unname(x$map),
					  stringsAsFactors = FALSE)
	coefs <- merge(coefs, map)
	coefs <- arrange(coefs, cM, strain)
	
	return(coefs)
	
}

#' Make a pretty-printing label for a locus with chr and pos
locus.label <- function(chr, pos, uu = NULL, ...) {
	
	## sniff units (Mb or cM)
	if (all(pos < 125) || uu == "cM") {
		pos <- sprintf("%.1f", pos)
		uu <- "cM"
	} else {
		pos <- sprintf("%.2f", pos/1e6)
		uu <- "Mb"
	}
	
	if (!any(grepl("^chr", chr))) {
		chr <- paste0("chr", chr)
	}
	
	paste0(chr, ": ", pos, " ", uu)
	
}

#' Perform genome scan and estimate allele effects at QTL peak(s).
scan.chrom <- function(G, phe, snps, chrom, cores = 8, ...) {
	
	this.geno <- subset(G, chr = chrom)
	lods <- scan1(this.geno, pheno, cores = cores)
	bi <- bayes_int(lods, chrom)
	map <- sort(lods$map[[1]])
	ii <- findInterval(bi, map)
	mk <- names(map)[ii]
	
	coefs <- scan1blup(grab.marker(this.geno, mk[2]), pheno, cores = cores)
	
	#df <- data.frame(marker = rownames(coefs$coef),
	#				 lod = lods$lod[ rownames(coefs$coef), ],
	#				 as.data.frame(coefs$coef))
	#df <- merge(df, snps)
	#df <- reshape2::melt(df, id.vars = c("chr","marker","cM","pos"))
	#df$what <- ifelse(df$variable == "lod", "LOD score", "coefficients")
	#df$variable <- factor(df$variable, LETTERS[1:8])
	
	return(list(lods = lods, coefs = coefs, best = mk[2], left = mk[1], right = mk[3]))
	
}

#' Tidy up the result of qtl2scan::scan1() for plotting LOD curves
tidy.scan1 <- function(x, ...) {
	
	## first unroll the marker map
	maps <- plyr::ldply(x$map, function(m) data.frame(marker = names(m), cM = unname(m), stringsAsFactors = FALSE), .id = "chr")
	rownames(maps) <- maps$marker
	
	## now get LOD scores
	lods <- reshape2::melt(x$lod)
	colnames(lods) <- c("marker","trait","lod")
	lods <- merge(lods, maps)
	lods <- plyr::arrange(lods, trait, chr, cM)
	
	return(lods)
	
}
