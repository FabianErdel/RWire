#' This function returns bootstrapped correlation coefficients
#'
#' @param am Accessibility matrix
#' @param replicates Number of replicates for the bootstrap
#' @return Table with correlation coefficients
#' @export

getCorRanges<-function(am, replicates = 100) {
  # check input arguments
  if(class(am)!="AccMatrix") stop("am must be an AccMatrix object")

  # determine number of ROIs and number of cells
  nrois <- dim(am@accmat)[1]
  ncells <- dim(am@accmat)[2]

  # calculate correlations for resampled accessibility matrices
  cortbl <- array(0, dim = c(nrois, nrois, 1+replicates))
  cortbl[, , 1] <- cor(t(am@accmat))
  res <- as.data.frame(cortbl[, , 1])

  for(i in 1:replicates) {
    indices <- sample.int(ncells, ncells, replace = TRUE)
    cortbl[, , 1+i] <- cor(t(am@accmat[,indices]))
  }

  # replace NA values that arise from empty rows
  cortbl[is.na(cortbl)] <- 2

  # sort correlation coefficients across samples (third index)
  cortbl <- apply(cortbl, c(1,2), sort)
  cortbl <- aperm(cortbl, c(2,3,1))

  # put back NA values
  cortbl[cortbl==2] <- NA

  # return correlation coefficients
  return(list(am@coord[,1], am@coord[,2], am@coord[,3], cortbl, res))
}
