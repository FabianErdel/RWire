#' This function makes cross correlation matrices from different sets of regions within accessibility matrices
#'
#' @param accmat1 Accessibility matrix 1 (data frame)
#' @param accmat2 Accessibility matrix 2 (data frame)
#' @param nmax Number of cells used for the analysis. When set to 0, all cells will be considered.
#' @return Cross correlation matrix
#' @export

makeCrossCorMatrix<-function(accmat1, accmat2, nmax = 0) {
  # check input arguments
  if(!is.data.frame(accmat1)) stop("accmat1 must be a data frame")
  if(!is.data.frame(accmat2)) stop("accmat2 must be a data frame")
  if(!is.numeric(nmax)) stop("nmax must be a number")

  # determine number of ROIs
  nrois1 <- dim(accmat1)[1]
  nrois2 <- dim(accmat2)[1]

  # determine number of cells
  ncells <- dim(accmat1)[2]-3

  # reduce ncells if applicable
  if((nmax > 0) & (ncells > nmax)) {ncells = nmax}

  # convert accessibility matrix to numeric matrix object
  am1 <- t(as.matrix(accmat1[1:nrois1, 4:(3+ncells)]))
  am2 <- t(as.matrix(accmat2[1:nrois2, 4:(3+ncells)]))

  # make correlation matrix
  cormat <- as.data.frame(cor(am1,am2))

  # add annotation columns
  cormat <- cbind(accmat1[,1], accmat1[,2], accmat1[,3], cormat)

  # add annotation rows
  nm <- as.data.frame(matrix(rep(c("chr", "start", "end"), 3), ncol=3, byrow=T), stringsAsFactors=F)
  nm <- as.data.frame(t(rbind(setNames(nm, names(accmat2[,1:3])), accmat2[,1:3])), stringsAsFactors=F)
  cormat <- rbind(nm, setNames(cormat, names(nm)))

  # return correlation matrix
  return(cormat)
}
