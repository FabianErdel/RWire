#' This function makes correlation matrices from accessibility matrices
#'
#' @param accmat Accessibility matrix (data frame)
#' @param nmax Number of cells used for the analysis. When set to 0, all cells will be considered.
#' @return Correlation matrix
#' @export

makeCorMatrix<-function(accmat, nmax = 0) {
  # check input arguments
  if(!is.data.frame(accmat)) stop("accmat must be a data frame")
  if(!is.numeric(nmax)) stop("nmax must be a number")

  # determine number of ROIs
  nrois <- dim(accmat)[1]

  # determine number of BED files in path
  ncells <- dim(accmat)[2]-3

  # reduce ncells if applicable
  if((nmax > 0) & (ncells > nmax)) {ncells = nmax}

  # convert accessibility matrix to numeric matrix object
  am <- t(as.matrix(accmat[1:nrois, 4:(3+ncells)]))

  # make correlation matrix
  cormat <- as.data.frame(cor(am))

  # add annotation columns
  cormat <- cbind(accmat[,1], accmat[,2], accmat[,3], cormat)

  # return correlation matrix
  return(cormat)
}
