#' This function returns p-values assessing the difference between two sets of correlation coefficients based on their bootstrapped confidence intervals
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

  # calculate confidence intervals for difference of correlation coefficients
  lower <- list()
  upper <- list()

  for(z in 1:((replicates-1)/2)) {
    lower[[z]] <- corrng1[[5]]-corrng2[[5]] - sqrt((corrng1[[5]] - corrng1[[4]][,,z])^2 + (corrng2[[4]][,,replicates+1-z] - corrng2[[5]])^2)
    upper[[z]] <- corrng1[[5]]-corrng2[[5]] + sqrt((corrng1[[4]][,,replicates+1-z] - corrng1[[5]])^2 + (corrng2[[5]] - corrng2[[4]][,,z])^2)
  }

  # convert lists into arrays
  lower <- array(as.numeric(unlist(lower)), dim=c(nrois, nrois, replicates))
  upper <- array(as.numeric(unlist(upper)), dim=c(nrois, nrois, replicates))

  # find highest confidence level, for which zero is not within the confidence interval
  pval <- matrix(1, nrow=nrois, ncol=nrois)

  for(x in 1:nrois) {
    for(y in 1:nrois) {
      if(sign(lower[x,y,1])==sign(upper[x,y,1]) & corrng1[[5]][x,y] != corrng2[[5]][x,y]) {
        pval[x,y] <- 2/(replicates-1)
      }
    }
  }

  for(x in 1:nrois) {
    for(y in 1:nrois) {
      for(z in 1:(replicates-1)) {
        if(sign(lower[x,y,z])!=sign(upper[x,y,z]) & sign(lower[x,y,z+1])==sign(upper[x,y,z+1])) {
          pval[x,y] <- 2*z/(replicates-1)
          break
        }
      }
    }
  }

  # return matrix with p-values
  return(list(corrng1[[1]],corrng1[[2]],corrng1[[3]],pval))
}
