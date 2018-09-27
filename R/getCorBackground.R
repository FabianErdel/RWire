#' This function returns the background correlation obtained with permuted entries
#'
#' @param am1 Accessibility matrix 1
#' @param am2 Accessibility matrix 2
#' @return Matrix with background correlations (quantiles)
#' @export

getCorBackground<-function(am1, am2) {
  # check input arguments
  if(class(am1)!="AccMatrix") stop("am1 must be an AccMatrix object")
  if(class(am2)!="AccMatrix") stop("am2 must be an AccMatrix object")
  if(dim(am1@accmat)[2]!=dim(am2@accmat)[2]) stop("dimensions are incompatible")

  # determine number of ROIs and number of cells
  nrois1 <- dim(am1@accmat)[1]
  nrois2 <- dim(am2@accmat)[1]
  ncells <- dim(am1@accmat)[2]

  # initialize data frame for results
  bg_cis<- matrix(0, nrow=2, ncol=100)
  bg_trans<- matrix(0, nrow=2, ncol=100)

  ## correlation matrix for resampled rows
  ## (total numbe of integrations for each genomic element kept constant)
  scram1 <- new("AccMatrix", coord = am1@coord, accmat = as.data.frame(t(apply(am1@accmat, 1, function(x) sample(x)))))
  cm <- makeCrossCorMatrix(scram1, am2)

  # get coordinates of the first element on each chromosome
  chr1 <- c(1, which((cm@coord1[1:(dim(cm@coord1)[1]-1),1] == cm@coord1[2:dim(cm@coord1)[1],1]) == F)+1, dim(cm@coord1)[1]+1)
  chr2 <- c(1, which((cm@coord2[1:(dim(cm@coord2)[2]-1),1] == cm@coord2[2:dim(cm@coord2)[2],1]) == F)+1, dim(cm@coord2)[2]+1)

  # aggregate cis and trans correlations
  cis <- vector("list", length(chr1)-1)
  trans <- vector("list", length(chr1)-1)
  for(i in 1:(length(chr1)-1)) {
    chrname <- cm@cormat[chr1[i],1]
    chrindex <- which(cm@cormat[1,chr2[1:(length(chr2)-1)]]==chrname)
    cis[[i]] <- as.matrix(cm@cormat[chr1[i]:(chr1[i+1]-1), chr2[chrindex]:(chr2[chrindex+1]-1)])
    trans[[i]] <- as.matrix(cm@cormat[chr1[i]:(chr1[i+1]-1), -c(1:3, chr2[chrindex]:(chr2[chrindex+1]-1))])
  }

  # determine quantiles
  bg_cis[1,] <- quantile(do.call(c,cis[!is.na(cis)]), seq(0.01,1,by=0.01), na.rm=T)
  bg_trans[1,] <- quantile(do.call(c,trans[!is.na(trans)]), seq(0.01,1,by=0.01), na.rm=T)

  ## correlation for resampled columns
  ## (total number of integrations for each cell kept constant)
  scram1 <- new("AccMatrix", coord = am1@coord, accmat = as.data.frame(apply(am1@accmat, 2, function(x) sample(x))))
  cm <- makeCrossCorMatrix(scram1, am2)

  # get coordinates of the first element on each chromosome
  chr1 <- c(1, which((cm@cormat[1:(dim(cm@cormat)[1]-1),1] == cm@cormat[2:dim(cm@cormat)[1],1]) == F)+1, dim(cm@cormat)[1]+1)
  chr2 <- c(1, which((cm@cormat[1,1:(dim(cm@cormat)[2]-1)] == cm@cormat[1,2:dim(cm@cormat)[2]]) == F)+1, dim(cm@cormat)[2]+1)

  # aggregate cis and trans correlations
  cis <- vector("list", length(chr1)-1)
  trans <- vector("list", length(chr1)-1)
  for(i in 1:(length(chr1)-1)) {
    chrname <- cormat[chr1[i],1]
    chrindex <- which(cm@cormat[1,chr2[1:(length(chr2)-1)]]==chrname)
    cis[[i]] <- as.matrix(cm@cormat[chr1[i]:(chr1[i+1]-1), chr2[chrindex]:(chr2[chrindex+1]-1)])
    trans[[i]] <- as.matrix(cm@cormat[chr1[i]:(chr1[i+1]-1), -c(1:3, chr2[chrindex]:(chr2[chrindex+1]-1))])
  }

  # determine quantiles
  bg_cis[2,] <- quantile(do.call(c,cis[!is.na(cis)]), seq(0.01,1,by=0.01), na.rm=T)
  bg_trans[2,] <- quantile(do.call(c,trans[!is.na(trans)]), seq(0.01,1,by=0.01), na.rm=T)

  # return correlation coefficients
  return(list(bg_cis, bg_trans))
}
