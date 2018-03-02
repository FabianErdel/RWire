#' Wrapper for making and storing scaled corrleation matrices within 1MB around a defined region
#'
#' @param accmat Accessibility matrix
#' @param anno_start Start of region of interest
#' @param anno_end End of region of interest
#' @param cutoff Use only regions that are visible in at least this number of cells. 
#' @param path Storage path for scaled correlation matrix
#' @export

writeMatrix<-function(accmat, anno_start, anno_end, cutoff, path) {
  # check input arguments
  if(!is.data.frame(accmat)) stop("accmat must be a data frame")
  if(!is.numeric(anno_start)) stop("anno_start must be a number")
  if(!is.numeric(anno_end)) stop("anno_end must be a number")
  
  # define displayed region
  start <- (anno_start+anno_end)/2 - 500000
  end <- (anno_start+anno_end)/2 + 500000
  
  # remove rows that contain only zero entries
  accmat <- accmat[rowSums(accmat[,4:dim(accmat)[2]])>0, ]
  
  # remove regions below visibility cutoff
  accmat <- accmat[rowSums(accmat[,4:dim(accmat)[2]] > 0)>cutoff*(dim(accmat)[2]-3)/100, ]
  
  # crop accessibility matrix to genomic region of interest
  amcrop <- accmat[accmat[,2]>start & accmat[,3]<end, ]
  
  # get number of accessible regions for each matrix
  nrois <- dim(amcrop)[1]-3
  
  # make correlation matrices
  cormat <- makeCorMatrix(amcrop)
  
  # make correlation matrices scaled according to genomic positions
  pngpath <- paste0(path, "_scaled_correlation_matrix.png")
  
  # determine image size; each pixel will corresponds to a genomic range of (end-start)/size
  size <- 200
  
  # make scaled matrices
  scm <- scaleCorMatrix(cormat, size)
  
  # delete left half of the matrix
  for(i in 1:size) {scm[i, 1:(i-1)] <- 0}
  
  # add annotation
  scale <- size/(end-start)
  for(i in ((anno_start-start)*scale):((anno_end-start)*scale)) {scm[i+2, i-2] <- -1}
  
  # save matrix
  saveImage(scm, pngpath, size, size, getCorLUT(), -1, 1)
}