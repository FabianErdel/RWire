#' Wrapper for making and storing scaled corrleation matrices within 1MB around a defined region
#'
#' @param am Accessibility matrix
#' @param chr Chromosome for region of interest
#' @param anno_start Start of region of interest
#' @param anno_end End of region of interest
#' @param vcutoff Use only regions that are visible in at least this percentage of cells
#' @param path Storage path for scaled correlation matrix
#' @export

writeMatrix<-function(am, chr, anno_start, anno_end, vcutoff, size, path) {
  # check input arguments
  if(class(am)!="AccMatrix") stop("am must be an AccMatrix object")
  if(!is.numeric(chr)) stop("chr must be a number")
  if(!is.numeric(anno_start)) stop("anno_start must be a number")
  if(!is.numeric(anno_end)) stop("anno_end must be a number")

  # define displayed region
  start <- (anno_start+anno_end)/2 - 500000
  end <- (anno_start+anno_end)/2 + 500000

  # remove rows that contain only zero entries
  am <- removeEmptyRows(am)

  # remove regions below visibility cutoff
  am <- removeRowsBelowCutoff(am, vcutoff)

  # crop accessibility matrix to genomic region of interest
  amcrop <- cropAccMatrix(am, chr, anno_start, anno_end)

  # get number of accessible regions for each matrix
  nrois <- dim(amcrop@accmat)[1]

  # make correlation matrices
  cormat <- makeCorMatrix(amcrop)

  # make correlation matrices scaled according to genomic positions
  pngpath <- paste0(path, "_scaled_correlation_matrix.png")

  # make scaled matrices; each pixel will corresponds to a genomic range of (end-start)/size
  scm <- scaleCorMatrix(cormat, chr, size)

  # delete left half of the matrix
  for(i in 1:size) {scm[i, 1:(i-1)] <- 0}

  # add annotation
  scale <- size/(end-start)
  for(i in ((anno_start-start)*scale):((anno_end-start)*scale)) {scm[i+2, i-2] <- -1}

  # save matrix
  saveImage(scm, pngpath, size, size, getCorLUT(), -1, 1)
}
