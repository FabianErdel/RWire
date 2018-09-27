#' This function makes a correlation matrix from two accessibility matrices
#'
#' @param am1 Accessibility matrix 1
#' @param am2 Accessibility matrix 2
#' @param nmax Number of cells used for the analysis. When set to 0, all cells will be considered.
#' @return Correlation matrix
#' @export

makeCorMatrix<-function(am1, am2 = am1, nmax = 0) {
  # check input arguments
  if(class(am1)!="AccMatrix") stop("am1 must be an AccMatrix object")
  if(class(am2)!="AccMatrix") stop("am2 must be an AccMatrix object")
  if(!is.numeric(nmax)) stop("nmax must be a number")
  if(dim(am1@accmat)[2]!=dim(am2@accmat)[2]) stop("dimensions are incompatible")

  # determine number of ROIs
  nrois1 <- dim(am1@accmat)[1]
  nrois2 <- dim(am2@accmat)[1]

  # determine number of cells
  ncells <- dim(am1@accmat)[2]

  # reduce ncells if applicable
  if((nmax > 0) & (ncells > nmax)) {ncells = nmax}

  # make correlation matrix
  cm <- new("CorMatrix", coord1 = am1@coord, coord2 = am2@coord, cormat = as.data.frame(cor(t(am1@accmat), t(am2@accmat))))

  # return correlation matrix
  return(cm)
}
