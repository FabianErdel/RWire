#' This function returns p-values assessing the difference between two bootstrap distributions of correlation coefficients based on the Kolmogorov-Smirnov test
#'
#' @param corrng1 Table 1 with annotations and correlation coefficients
#' @param corrng2 Table 2 with annotations and correlation coefficients
#' @return Matrix with p-values
#' @export

compareCorDistributions<-function(corrng1, corrng2) {
  # check input arguments
  if(!is.list(corrng1)) stop("corrng1 must be a list")
  if(!is.list(corrng2)) stop("corrng2 must be a list")

  # retrieve array dimensions
  nrois <- dim(corrng1[[4]])[1]

  # perform KS test
  pval <- matrix(1, nrow=nrois, ncol=nrois)

  for(x in 1:nrois) {
    for(y in 1:nrois) {
      pval[x,y] <- ks.test(corrng1[[4]][x,y,], corrng2[[4]][x,y,])$p
    }
  }

  # return matrix with (approximated) p-values
  return(list(corrng1[[1]],corrng1[[2]],corrng1[[3]],pval))
}
