#' This function returns p-values assessing the difference between two sets of correlation coefficients based on the Kolmogorov-Smirnov test
#'
#' @param corrng1 Table 1 with annotations and correlation coefficients
#' @param corrng2 Table 2 with annotations and correlation coefficients
#' @return Matrix with p-values
#' @export

compareCorrelations<-function(corrng1, corrng2) {
  # check input arguments
  if(!is.list(corrng1)) stop("corrng1 must be a list")
  if(!is.list(corrng2)) stop("corrng2 must be a list")

  # retrieve array dimensions
  nrois <- dim(corrng1[[4]])[1]

  # perform KS test
  funks <- function(x,y) {
    return(ks.test(corrng1[[4]][x,y,], corrng2[[4]][x,y,])$p)
  }

  pval <- sapply(1:nrois, function(x) mapply(funks, x, 1:nrois))

  # return matrix with (approximated) p-values
  return(pval)
}
