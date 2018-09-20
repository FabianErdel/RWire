#' This function returns the background correlation obtained with permuted entries
#'
#' @param accmat1 Accessibility matrix 1
#' @param accmat2 Accessibility matrix 2
#' @return Matrix with background correlations (quantiles)
#' @export

getCorBackground<-function(accmat1, accmat2) {
  # check input arguments
  if(!is.data.frame(accmat1)) stop("accmat1 must be a data frame")
  if(!is.data.frame(accmat2)) stop("accmat2 must be a data frame")
  if(dim(accmat1)[2]!=dim(accmat2)[2]) stop("dimensions are incompatible")

  # determine number of ROIs and number of cells
  nrois1 <- dim(accmat1)[1]
  nrois2 <- dim(accmat2)[2]
  ncells <- dim(accmat1)[2]-3

  # convert data frames into matrices (omitting annotation columns)
  m1 <- data.matrix(accmat1[,4:dim(accmat1)[2]])
  m2 <- data.matrix(accmat2[,4:dim(accmat2)[2]])

  # initialize data frame for results
  bg_cis<- matrix(0, nrow=2, ncol=100)
  bg_trans<- matrix(0, nrow=2, ncol=100)

  ## correlation for resampled rows
  ## (total numbe of integrations for each genomic element kept constant)
  cormat <- makeCrossCorMatrix(t(apply(m1, 1, function(x) sample(x))), m2)

  # get coordinates of the first element on each chromosome
  chr1 <- c(4, which((cormat[4:(dim(cormat)[1]-1),1] == cormat[5:dim(cormat)[1],1]) == F)+4, dim(cormat)[1]+1)
  chr2 <- c(4, which((cormat[1,4:(dim(cormat)[2]-1)] == cormat[1,5:dim(cormat)[2]]) == F)+4, dim(cormat)[2]+1)

  # aggregate cis and trans correlations
  cis <- vector("list", length(chr1)-1)
  trans <- vector("list", length(chr1)-1)
  for(i in 1:(length(chr1)-1)) {
    chrname <- cormat[chr1[i],1]
    chrindex <- which(cormat[1,chr2[1:(length(chr2)-1)]]==chrname)
    cis[[i]] <- apply(as.matrix(cormat[chr1[i]:(chr1[i+1]-1), chr2[chrindex]:(chr2[chrindex+1]-1)]), c(1,2), as.numeric)
    trans[[i]] <- apply(as.matrix(cormat[chr1[i]:(chr1[i+1]-1), -c(1:3, chr2[chrindex]:(chr2[chrindex+1]-1))]), c(1,2), as.numeric)
  }

  # determine quantiles
  bg_cis[1,] <- quantile(do.call(c,cis[!is.na(cis)]), seq(0.01,1,by=0.01))
  bg_trans[1,] <- quantile(do.call(c,trans[!is.na(trans)]), seq(0.01,1,by=0.01))

  ## correlation for resampled columns
  ## (total number of integrations for each cell kept constant)
  cormat <- makeCrossCorMatrix(apply(m1, 2, function(x) sample(x)), m2)

  # get coordinates of the first element on each chromosome
  chr1 <- c(4, which((cormat[4:(dim(cormat)[1]-1),1] == cormat[5:dim(cormat)[1],1]) == F)+4, dim(cormat)[1]+1)
  chr2 <- c(4, which((cormat[1,4:(dim(cormat)[2]-1)] == cormat[1,5:dim(cormat)[2]]) == F)+4, dim(cormat)[2]+1)

  # aggregate cis and trans correlations
  cis <- vector("list", length(chr1)-1)
  trans <- vector("list", length(chr1)-1)
  for(i in 1:(length(chr1)-1)) {
    chrname <- cormat[chr1[i],1]
    chrindex <- which(cormat[1,chr2[1:(length(chr2)-1)]]==chrname)
    cis[[i]] <- apply(as.matrix(cormat[chr1[i]:(chr1[i+1]-1), chr2[chrindex]:(chr2[chrindex+1]-1)]), c(1,2), as.numeric)
    trans[[i]] <- apply(as.matrix(cormat[chr1[i]:(chr1[i+1]-1), -c(1:3, chr2[chrindex]:(chr2[chrindex+1]-1))]), c(1,2), as.numeric)
  }

  # determine quantiles
  bg_cis[2,] <- quantile(do.call(c,cis[!is.na(cis)]), seq(0.01,1,by=0.01))
  bg_trans[2,] <- quantile(do.call(c,trans[!is.na(trans)]), seq(0.01,1,by=0.01))

  # return correlation coefficients
  return(list(bg_cis, bg_trans))
}
