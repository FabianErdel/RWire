#' This function counts reads in a GRanges object that fall into a list of ROIs
#'
#' @param data GRanges object containing reads
#' @param rois GRanges object containing ROIs
#' @param mino Minimum overlap (see package 'GenomicRanges')
#' @param ignore_strand When set to TRUE, the strand information is ignored in the overlap calculations (see package 'GenomicRanges')
#' @param chr Chromosome number. When set to 0, all chromosomes will be considered.
#' @return GRanges object containing read counts in ROIs
#' @export

countReads<-function(data, rois, mino = 15, ignore_strand = TRUE, chr = 0) {
  # check input arguments
  if(!is.numeric(chr)) stop("chr must be a number")
  
  # filter rois for selected chromosome
  if(chr > 0) {rois <- rois[seqnames(rois) == paste0('chr',chr)]} 
    
  # count reads in ROIs
  counts <- GenomicRanges::countOverlaps(rois, data, minoverlap=mino, type="any", ignore.strand=ignore_strand)
  
  # return counts
  return(counts)
}