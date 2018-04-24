#' This function returns p-values assessing the difference between two sets of correlation coefficients based on their bootstrap distributions
#'
#' @param corrng1 Table 1 with annotations and correlation coefficients
#' @param corrng2 Table 2 with annotations and correlation coefficients
#' @return Matrix with p-values
#' @export

compareCorCoeffs<-function(corrng1, corrng2) {
  # check input arguments
  if(!is.list(corrng1)) stop("corrng1 must be a list")
  if(!is.list(corrng2)) stop("corrng2 must be a list")

  # retrieve array dimensions
  nrois <- dim(corrng1[[4]])[1]
  replicates <- dim(corrng1[[4]])[3]

  # retrieve p-values from bootstrap distribution of the difference of correlation coefficients
  pval <- matrix(1, nrow=nrois, ncol=nrois)

  for(x in 1:nrois) {
    for(y in 1:nrois) {
      # generate bootstrap distribution for difference of correlation coefficients
      cmb <- expand.grid(corrng1[[4]][x,y,], corrng2[[4]][x,y,])
      cmb <- cmb[is.finite(cmb[,1]) & is.finite(cmb[,2]),]

      # caclulate p-value (if valid bootstrap samples are available)
      if(dim(cmb)[1]>0) {
        if(corrng1[[5]][x,y]-corrng2[[5]][x,y]>0) {pval[x,y] <- (sum(cmb[,1]-cmb[,2]<0))/dim(cmb)[1]}
        if(corrng1[[5]][x,y]-corrng2[[5]][x,y]<0) {pval[x,y] <- (sum(cmb[,1]-cmb[,2]>0))/dim(cmb)[1]}
      }
    }
  }

  # set the diagonal to unity
  for(i in 1:nrois) {
    pval[i,i] <- 1
  }

  # return matrix with p-values
  return(list(corrng1[[1]],corrng1[[2]],corrng1[[3]],pval))
}
