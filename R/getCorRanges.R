#' This function returns correlation coefficients and their confidence intervals
#'
#' @param accmat Accessibility matrix
#' @param replicates Number of replicates for the bootstrap
#' @param confidence Confidence level (between 0 and 1)
#' @return Stack of matrices containing correlation coefficients (first layer) along with minimum, maximum and median values from bootstrapping (second, third and fourth layer)
#' @export

getCorRanges<-function(accmat, replicates = 100, confidence = 0.95) {
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
  res1 <- as.data.frame(cortbl[, , 1])
  
  for(i in 1:replicates) {
    indices <- sample.int(ncells, ncells, replace = TRUE)
    cortbl[, , 1+i] <- cor(t(m[,indices]))
  }
  
  # remove NA values that arise from empty rows
  cortbl[is.na(cortbl)] <- 0

  # sort correlation coefficients across samples (third index)
  cortbl <- apply(cortbl, c(1,2), sort)
  cortbl <- aperm(cortbl, c(2,3,1))

  # retrieve start/end of confidence intervals and median correlation coefficient
  tail <- round(0.5*(1-confidence)*replicates)
  start_index <- tail
  end_index <- replicates - tail + 2

  res2 <- as.data.frame(cortbl[, , start_index])
  res3 <- as.data.frame(cortbl[, , end_index])
  res4 <- as.data.frame(cortbl[, , 1+replicates/2])

  # add annotation columns
  res1 <- cbind(accmat[,1], accmat[,2], accmat[,3], res1) # correlation coefficient without bootstrapping
  res2 <- cbind(accmat[,1], accmat[,2], accmat[,3], res2) # minimum correlation coefficient from bootstrapping
  res3 <- cbind(accmat[,1], accmat[,2], accmat[,3], res3) # maximum correlation coefficient from bootstrapping
  res4 <- cbind(accmat[,1], accmat[,2], accmat[,3], res4) # median correlation coefficient from bootstrapping
  
  # return correlation coefficients and start/end positions of confidence intervals
  return(list(res1,res2,res3,res4))
}