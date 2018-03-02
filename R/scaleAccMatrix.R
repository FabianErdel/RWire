#' This function scales accessibility matrices to genomic coordinates
#'
#' @param accmat Accessibility matrix
#' @param size Number of genomic regions to be used
#' @param start Minimum genomic coordinate
#' @param end Maximum genomic coordinate
#' @return Scaled accessibility matrix
#' @export

scaleAccMatrix<-function(accmat, size = 1000, start = 0, end = 249250621) {
  # check input arguments
  if(!is.data.frame(accmat)) stop("accmat must be a data frame")
  if(!is.numeric(size)) stop("size must be a number")
  if(!is.numeric(start)) stop("start must be a number")
  if(!is.numeric(end)) stop("end must be a number")

  # determine number of ROIs
  nrois <- dim(accmat)[1]
  ncells <- dim(accmat)[2]-3

  # make numeric accessibility matrix containing regions between start and end
  sel_indices <- accmat[,2]>start & accmat[,3]<end
  m <- data.matrix(accmat[sel_indices, 2:dim(accmat)[2]])

  # determine scaling factor (kb/pixel)
  # scale <- (max(m[,1])-min(m[,1]))/size

  # make genomic bins
  counts <- as.numeric(table(cut(m[,1], size)))
  indices <- c(0, cumsum(counts))
  
  # make binned accessibility matrix
  bm <- matrix(0, ncol = ncells, nrow = size)

  for(i in 1:size) {
    if(counts[i] > 0) {
      bm[i, ] <- colSums(m[(indices[i]+1):indices[i+1], ])[3:(2+ncells)]
    }
  }

  # return binned accessibility matrix
  return(bm)
}