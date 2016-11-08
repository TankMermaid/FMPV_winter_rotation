# wgs.R
# some utility functions related to whole-genome sequencing

#' Read output from my \code{coverage.py} script
#' 
#' @param path name of file to read
#' @param scrub if \code{TRUE}, remove trailing suffixes (\code{*.sorted}, \code{*.merged} etc) from individual IDs
#' @param as.ranges if \code{TRUE}, convert to a \code{GRanges} object
#' 
#' @details Expected input is a BED-like file with columns *chr*, *start*, *end*, *iid*, *raw read count*,
#' 	*normalized read count* and optionally *count of MQ0 reads*, *proportion of total reads which are MQ0*.
#' 
#' @return a dataframe
read.coverage <- function(path, scrub = TRUE, as.ranges = FALSE, has.MQ0 = FALSE, ...) {
	
	counts <- tryCatch({
		readr::read_tsv(path, col_names = FALSE, comment = "#")
	},
	error = function(e) {
		read.table(path, header = FALSE, comment.char = "#")
	})
	colnames(counts)[1:6] <- c("chr","start","end","iid","nreads","depth")
	if (has.MQ0)
		colnames(counts)[7:8] <- c("MQ0","MQ0.prop")
	
	if (scrub)
		counts$iid <- gsub("\\.\\w+$", "", counts$iid)
	
	if (as.ranges)
		counts <- GenomicRanges::makeGRangesFromDataFrame(counts, starts.in.df.are.0based = TRUE,
														  keep.extra.columns = TRUE)
	
	attr(counts, "filename") <- normalizePath(path)
	return(counts)
	
}

read.bed <- function(path, header = FALSE, as.ranges = FALSE, seqinfo = NULL, ...) {
	
	bed <- tryCatch({
		readr::read_tsv(path, col_names = header, comment = "#")
	},
	error = function(e) {
		read.table(path, header = header, sep = "\t")
	})
	
	cols <- c("chr","start","end","name","score","strand")
	colnames(bed)[ 1:min(ncol(bed), length(cols)) ] <- cols[ 1:ncol(bed) ]
	
	if (as.ranges) {
		bed <- GenomicRanges::makeGRangesFromDataFrame(bed, starts.in.df.are.0based = TRUE,
													   keep.extra.columns = TRUE, seqinfo = seqinfo)
	} else {
		if (!is.null(seqinfo))
			bed$chr <- factor(bed$chr, seqlevels(seqinfo))
	}
	
	attr(bed, "filename") <- normalizePath(path)
	return(bed)
	
}

write.bed <- function(gr, outfile, extra = FALSE, zero.based = TRUE, ...) {
	
	df <- as.data.frame(gr)
	ii <- which(colnames(df) == "seqnames")
	if (length(ii)) {
		colnames(df)[ii[1]] <- "chr"
	}
	
	if (inherits(gr, "GRanges") || !zero.based)
		df$start <- pmax(0, df$start-1)
	
	towrite <- df[ ,c("chr","start","end") ]
	if (extra)
		towrite <- cbind(towrite, df[ ,setdiff(colnames(df), c("chr","start","end")) ])
	
	write.table( towrite, outfile, quote = FALSE, sep = "\t",
				 row.names = FALSE, col.names = FALSE )
	
}

read.picard.metrics <- function(fname, ...) {
	
	ll <- readLines(fname)
	bstart <- grep("^## [A-Z]+", ll)
	blanks <- which(nchar(ll) == 0)
	bpos <- findInterval(blanks, bstart)
	bend <- blanks[ bpos > 0 ]-1
	blocks <- cbind(start = bstart, end = bend)
	
	rez <- plyr::alply(blocks, 1, function(row) {
		streamer <- textConnection(ll[ row[1]:row[2] ])
		read.table(streamer, header = TRUE)
	})
	
	bnames <- ll[bstart]
	bnames <- gsub("^## ", "", bnames)
	names(rez) <- bnames
	return(rez)
	
}