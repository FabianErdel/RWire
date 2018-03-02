#' This function makes tiles for a particular chromosome of the hg19 genome
#'
#' @param chr Chromosome number
#' @param start Start position (if only tiles for a particular region should be made)
#' @param end End position (if only tiles for a particular region should be made) 
#' @param tile_size Tile size
#' @param ignore_strand When set to TRUE, the strand information is ignored in the overlap calculations (see package 'GenomicRanges')
#' @return GRanges object with tiles for the indicated chromosome
#' @export

makeTiles<-function(chr = 1, start = 0, end = 0, tile_size = 1000, ignore_strand = FALSE) {
  # check input arguments
  if(!is.numeric(tile_size)) stop("tile size must be a number")
  if(!is.numeric(start)) stop("start must be a number")
  if(!is.numeric(end)) stop("end must be a number")
  if(!is.numeric(chr)) stop("chr must be a number")
  
  # adjust end for each chromosome
  if(end == 0) {
    if(chr == 1) {end = 249250621}
    if(chr == 2) {end = 243199373}
    if(chr == 3) {end = 198022430}
    if(chr == 4) {end = 191154276}
    if(chr == 5) {end = 180915260}
    if(chr == 6) {end = 171115067}
    if(chr == 7) {end = 159138663}
    if(chr == 8) {end = 146364022}
    if(chr == 9) {end = 141213431}
    if(chr == 10) {end = 135534747}
    if(chr == 11) {end = 135006516}
    if(chr == 12) {end = 133851895}
    if(chr == 13) {end = 115169878}
    if(chr == 14) {end = 107349540}
    if(chr == 15) {end = 102531392}
    if(chr == 16) {end = 90354753}
    if(chr == 17) {end = 81195210}
    if(chr == 18) {end = 78077248}
    if(chr == 19) {end = 59128983}
    if(chr == 20) {end = 63025520}
    if(chr == 21) {end = 48129895}
    if(chr == 22) {end = 51304566}
    
    if(chr == 23) {end = 155270560}
    if(chr == 24) {end = 59373566}
  }
  
  # make data frame with tiles
  df <- data.frame(chr=paste0("chr",chr), start=seq(start, end-tile_size, tile_size) + 1, end=seq(start, end-tile_size, tile_size) + tile_size,
        strand="+")
  
  # assign column names for conversion to GRanges object
  colnames(df) <- c("chr", "start", "end", "strand")
  
  # convert data table to GRanges object
  tiles <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=FALSE, ignore.strand=ignore_strand, seqinfo=NULL,
                                    seqnames.field="chr", start.field="start", end.field="end", strand.field="strand", starts.in.df.are.0based=FALSE)
  
  # return tiles
  return(tiles)
}