#' This function returns bootstrapped correlation coefficients
#'
#' @param accmat Accessibility matrix
#' @param replicates Number of replicates for the bootstrap
#' @return Table with correlation coefficients
#' @export

getCorRanges<-function(accmat, replicates = 100) {
  # check input arguments
  if(!is.data.frame(accmat)) stop("accmat must be a data frame")

  # determine number of ROIs and number of cells
  nrois <- dim(accmat)[1]
  ncells <- dim(accmat)[2]-3

  # convert data frame into matrix (omitting annotation columns)
  m <- data.matrix(accmat[,4:dim(accmat)[2]])

  # calculate correlations for resampled accessibility matrices
  cortbl <- array(0, dim = c(nrois, nrois, 1+replicates))
  cortbl[, , 1] <- cor(t(m))
  res <- as.data.frame(cortbl[, , 1])

  for(i in 1:replicates) {
    indices <- sample.int(ncells, ncells, replace = TRUE)
    cortbl[, , 1+i] <- cor(t(m[,indices]))
  }

  # replace NA values that arise from empty rows
  cortbl[is.na(cortbl)] <- 2

  # sort correlation coefficients across samples (third index)
  cortbl <- apply(cortbl, c(1,2), sort)
  cortbl <- aperm(cortbl, c(2,3,1))

  # put back NA values
  cortbl[cortbl==2] <- NA

  # return correlation coefficients
  return(list(accmat[,1], accmat[,2], accmat[,3], cortbl, res))
}
