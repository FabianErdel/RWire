#' This function makes tiles for the complete hg19 genome
#'
#' @param tile_size Tile size
#' @return GRanges object with tiles for the complete genome
#' @export

makeGenomeTiles<-function(tile_size) {
  # check input arguments
  if(!is.numeric(tile_size)) stop("tile size must be a number")
  
  # make data frames with tiles
  df <- list()
  df[[1]] <- data.frame(chr="chr1", start=seq(0, 249250621, tile_size) + 1, end=seq(0, 249250621, tile_size) + tile_size, strand="+")
  df[[2]] <- data.frame(chr="chr2", start=seq(0, 243199373, tile_size) + 1, end=seq(0, 243199373, tile_size) + tile_size, strand="+")
  df[[3]] <- data.frame(chr="chr3", start=seq(0, 198022430, tile_size) + 1, end=seq(0, 198022430, tile_size) + tile_size, strand="+")
  df[[4]] <- data.frame(chr="chr4", start=seq(0, 191154276, tile_size) + 1, end=seq(0, 191154276, tile_size) + tile_size, strand="+")
  df[[5]] <- data.frame(chr="chr5", start=seq(0, 180915260, tile_size) + 1, end=seq(0, 180915260, tile_size) + tile_size, strand="+")
  df[[6]] <- data.frame(chr="chr6", start=seq(0, 171115067, tile_size) + 1, end=seq(0, 171115067, tile_size) + tile_size, strand="+")
  df[[7]] <- data.frame(chr="chr7", start=seq(0, 159138663, tile_size) + 1, end=seq(0, 159138663, tile_size) + tile_size, strand="+")
  df[[8]] <- data.frame(chr="chr8", start=seq(0, 146364022, tile_size) + 1, end=seq(0, 146364022, tile_size) + tile_size, strand="+")
  df[[9]] <- data.frame(chr="chr9", start=seq(0, 141213431, tile_size) + 1, end=seq(0, 141213431, tile_size) + tile_size, strand="+")
  df[[10]] <- data.frame(chr="chr10", start=seq(0, 135534747, tile_size) + 1, end=seq(0, 135534747, tile_size) + tile_size, strand="+")
  df[[11]] <- data.frame(chr="chr11", start=seq(0, 135006516, tile_size) + 1, end=seq(0, 135006516, tile_size) + tile_size, strand="+")
  df[[12]] <- data.frame(chr="chr12", start=seq(0, 133851895, tile_size) + 1, end=seq(0, 133851895, tile_size) + tile_size, strand="+")
  df[[13]] <- data.frame(chr="chr13", start=seq(0, 115169878, tile_size) + 1, end=seq(0, 115169878, tile_size) + tile_size, strand="+")
  df[[14]] <- data.frame(chr="chr14", start=seq(0, 107349540, tile_size) + 1, end=seq(0, 107349540, tile_size) + tile_size, strand="+")
  df[[15]] <- data.frame(chr="chr15", start=seq(0, 102531392, tile_size) + 1, end=seq(0, 102531392, tile_size) + tile_size, strand="+")
  df[[16]] <- data.frame(chr="chr16", start=seq(0, 90354753, tile_size) + 1, end=seq(0, 90354753, tile_size) + tile_size, strand="+")
  df[[17]] <- data.frame(chr="chr17", start=seq(0, 81195210, tile_size) + 1, end=seq(0, 81195210, tile_size) + tile_size, strand="+")
  df[[18]] <- data.frame(chr="chr18", start=seq(0, 78077248, tile_size) + 1, end=seq(0, 78077248, tile_size) + tile_size, strand="+")
  df[[19]] <- data.frame(chr="chr19", start=seq(0, 59128983, tile_size) + 1, end=seq(0, 59128983, tile_size) + tile_size, strand="+")
  df[[20]] <- data.frame(chr="chr20", start=seq(0, 63025520, tile_size) + 1, end=seq(0, 63025520, tile_size) + tile_size, strand="+")
  df[[21]] <- data.frame(chr="chr21", start=seq(0, 48129895, tile_size) + 1, end=seq(0, 48129895, tile_size) + tile_size, strand="+")
  df[[22]] <- data.frame(chr="chr22", start=seq(0, 51304566, tile_size) + 1, end=seq(0, 51304566, tile_size) + tile_size, strand="+")
  
  df[[23]] <- data.frame(chr="chrX", start=seq(0, 155270560, tile_size) + 1, end=seq(0, 155270560, tile_size) + tile_size, strand="+")
  df[[24]] <- data.frame(chr="chrY", start=seq(0, 59373566, tile_size) + 1, end=seq(0, 59373566, tile_size) + tile_size, strand="+")

  # combine in one list
  all <- rbindlist(df)
  
  # assign column names for conversion to GRanges object
  colnames(all) <- c("chr", "start", "end", "strand")
  
  # convert data table to GRanges object
  tiles <- GenomicRanges::makeGRangesFromDataFrame(all, keep.extra.columns=FALSE, ignore.strand=FALSE, seqinfo=NULL,
           seqnames.field="chr", start.field="start", end.field="end", strand.field="strand", starts.in.df.are.0based=FALSE)
  
  # return tiles
  return(tiles)
}