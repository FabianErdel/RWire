#' This function returns promoter-enhancer (P-E) pairs based on a correlation matrix
#'
#' @param cm Correlation matrix (rows: promoters, columns: enhancers)
#' @param cutoff Cutoff value required for P-E assignment
#' @param window Genomic window for P-E assignment
#' @return P-E pairs
#' @export

getPEPairs<-function(cm, cutoff, window = .Machine$double.xmax) {
  # check input arguments
  if(class(cm)!="CorMatrix") stop("cm must be a CorMatrix object")
  if(!is.numeric(cutoff)) stop("cutoff must be a number")
  if(!is.numeric(window)) stop("window must be a number")

  # create PEPairs object
  pep <- new("PEPairs", promoters = cm@coord1, enhancers = cm@coord2)
  pep@promoters[,dim(pep@promoters)[2]+1] <- rep(0, dim(pep@promoters)[1])
  pep@enhancers[,dim(pep@enhancers)[2]+1] <- rep(0, dim(pep@enhancers)[1])
  colnames(pep@promoters)[dim(pep@promoters)[2]] <- "assigned enhancer indices"
  colnames(pep@enhancers)[dim(pep@enhancers)[2]] <- "assigned promoter indices"

  # assign enhancers to each promoter
  for(i in 1:dim(cm@coord1)[1]) { # loop through promoters
    # select enhancers within genomic window
    #d <- abs((cm@coord1[i,2]+cm@coord1[i,3])/2-(cm@coord2[,2]+cm@coord2[,3])/2)
    elist <- which(cm@coord1[i,1]==cm@coord2[,1] & abs((cm@coord1[i,2]+cm@coord1[i,3])/2-(cm@coord2[,2]+cm@coord2[,3])/2)>0 & abs((cm@coord1[i,2]+cm@coord1[i,3])/2-(cm@coord2[,2]+cm@coord2[,3])/2)<=window)

    if(length(elist)>0) {
      # find enhancer above cutoff
      indices <- which(cm@cormat[i,elist] >= cutoff)
      pep@promoters$`assigned enhancer indices`[i] <- list(elist[indices])
    }
  }

  # assign promoters to each enhancer
  for(i in 1:dim(cm@coord2)[1]) { # loop through enhancers
    # select promoters within genomic window
    #d <- abs((cm@coord1[i,2]+cm@coord1[i,3])/2-(cm@coord2[,2]+cm@coord2[,3])/2)
    plist <- which(cm@coord1[i,1]==cm@coord2[,1] & abs((cm@coord1[i,2]+cm@coord1[i,3])/2-(cm@coord2[,2]+cm@coord2[,3])/2)>0 & abs((cm@coord1[i,2]+cm@coord1[i,3])/2-(cm@coord2[,2]+cm@coord2[,3])/2)<=window)

    if(length(plist)>0) {
      # find enhancer above cutoff
      indices <- which(cm@cormat[plist,i] >= cutoff)
      pep@enhancers$`assigned promoter indices`[i] <- list(plist[indices])
    }
  }

  # return correlation coefficients
  return(pep)
}
