#' This function reads a bed file and stores the reads in a GRanges object
#'
#' @param filename File name of bed file
#' @param chr Chromosome number
#' @param ignore_strand When set to TRUE, the strand information is ignored in the overlap calculations (see package 'GenomicRanges')
#' @return GRanges object with reads
#' @export

readBed<-function(filename, chr = 1, ignore_strand = FALSE) {
  # check input arguments
  if(!file.exists(filename)) stop("file not found")
  if(!is.numeric(chr)) stop("chr must be a number")

  # read data
  data <- data.table::fread(filename, data.table=FALSE, header=FALSE)

  # assign column names for conversion to GRanges object
  colnames(data) <- c("chr", "start", "end", "name", "score", "strand")

  # convert data table to GRanges object
  reads <- GenomicRanges::makeGRangesFromDataFrame(data, keep.extra.columns=FALSE, ignore.strand=ignore_strand, seqinfo=NULL,
       seqnames.field="chr", start.field="start", end.field="end", strand.field="strand", starts.in.df.are.0based=FALSE)

  # return data
  return(reads)
}
