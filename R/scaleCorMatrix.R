#' This function scales correlation matrices to genomic coordinates
#'
#' @param cormat Correlation matrix
#' @param size Number of genomic regions to be used
#' @param start Minimum genomic coordinate
#' @param end Maximum genomic coordinate
#' @return Scaled correlation matrix
#' @export

scaleCorMatrix<-function(cormat, size = 1000, start = 0, end = 249250621) {
  # check input arguments
  if(!is.data.frame(cormat)) stop("cormat must be a data frame")
  if(!is.numeric(size)) stop("size must be a number")
  if(!is.numeric(start)) stop("start must be a number")
  if(!is.numeric(end)) stop("end must be a number")

  # determine number of ROIs
  nrois <- dim(cormat)[1]

  # make numeric correlation matrix containing regions between start and end
  sel_indices <- cormat[,2]>start & cormat[,3]<end
  m <- data.matrix(cormat[sel_indices, c(FALSE, TRUE, TRUE, sel_indices)])

  # determine scaling factor (kb/pixel)
  # scale <- (max(m[,1])-min(m[,1]))/size

  # make genomic bins
  counts <- as.numeric(table(cut(m[,1], size)))
  indices <- c(0, cumsum(counts))

  # make binned correlation matrix
  bm <- matrix(0, ncol = size, nrow = size)
    
  for(i in 1:size) {
    for(j in 1:size) {
      if(counts[i]) {
        bm[i,j] <- mean(m[(indices[i]+1):indices[i+1], (2+indices[j]+1):(2+indices[j+1])])
      }
    }
  }
  
  # return binned correlation matrix
  return(bm)
}