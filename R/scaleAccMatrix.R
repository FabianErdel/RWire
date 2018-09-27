#' This function scales accessibility matrices to genomic coordinates
#'
#' @param am Accessibility matrix
#' @param chr Chromosome to be used
#' @param size Number of genomic regions to be used
#' @param start Minimum genomic coordinate
#' @param end Maximum genomic coordinate
#' @return Scaled accessibility matrix
#' @export

scaleAccMatrix<-function(am, chr, size = 1000, start = 0, end = 249250621) {
  # check input arguments
  if(class(am)!="AccMatrix") stop("am must be an AccMatrix object")
  if(!is.numeric(chr)) stop("chr must be a number")
  if(!is.numeric(size)) stop("size must be a number")
  if(!is.numeric(start)) stop("start must be a number")
  if(!is.numeric(end)) stop("end must be a number")

  # determine number of ROIs
  nrois <- dim(am@accmat)[1]
  ncells <- dim(am@accmat)[2]

  # crop accessibility matrix
  am <- cropAccMatrix(am, chr, start, end)

  # determine scaling factor (kb/pixel)
  # scale <- (max(am@coord[,2])-min(am@coord[,2]))/size

  # make genomic bins
  counts <- as.numeric(table(cut(am@coord[,2], size)))
  indices <- c(0, cumsum(counts))

  # make binned accessibility matrix
  bm <- matrix(0, ncol = ncells, nrow = size)

  for(i in 1:size) {
    if(counts[i] > 0) {
      bm[i, ] <- colSums(am@accmat[(indices[i]+1):indices[i+1], ])
    }
  }

  # return binned accessibility matrix
  return(bm)
}
