#' This function scales correlation matrices to genomic coordinates
#'
#' @param cm Correlation matrix
#' @param chr Number of chromosome
#' @param size Number of genomic regions to be used
#' @param start Minimum genomic coordinate
#' @param end Maximum genomic coordinate
#' @return Scaled correlation matrix
#' @export

scaleCorMatrix<-function(cm, chr, size = 1000, start = 0, end = 249250621) {
  # check input arguments
  if(class(cm)!="CorMatrix") stop("cm must be a CorMatrix object")
  if(!identical(cm@coord1,cm@coord2)) stop("cm must be an autocorrelation matrix")
  if(!is.numeric(chr)) stop("chr must be a number")
  if(!is.numeric(size)) stop("size must be a number")
  if(!is.numeric(start)) stop("start must be a number")
  if(!is.numeric(end)) stop("end must be a number")

  # determine number of ROIs
  nrois <- dim(cm@cormat)[1]

  # make numeric correlation matrix containing regions (on chromosome chr) between start and end
  sel_indices <- cm@coord1[,1]==paste0('chr',chr) & cm@coord1[,2]>start & cm@coord1[,3]<end
  m <- data.matrix(cm@cormat[sel_indices, sel_indices])

  # determine scaling factor (kb/pixel)
  # scale <- (max(m[,1])-min(m[,1]))/size

  # make genomic bins
  counts <- as.numeric(table(cut(m[,1], size)))
  indices <- c(0, cumsum(counts))

  # make binned correlation matrix
  bm <- matrix(0, ncol = size, nrow = size)

  for(i in 1:size) {
    for(j in 1:size) {
      if(counts[i] & counts[j]) {
        bm[i,j] <- mean(m[(indices[i]+1):indices[i+1], (2+indices[j]+1):(2+indices[j+1])])
      }
    }
  }

  # return binned correlation matrix
  return(bm)
}
