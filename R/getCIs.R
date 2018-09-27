#' This function returns confidence intervals
#'
#' @param corrng Table with annotations and correlation coefficients
#' @param confidence Confidence level (between 0 and 1)
#' @return Stack of CorMatrices containing correlation coefficients (first layer) along with minimum, maximum and median values from bootstrapping (second, third and fourth layer)
#' @export

getCIs<-function(corrng, confidence = 0.95) {
  # check input arguments
  if(!is.list(corrng)) stop("corrng must be a list")

  # treat NA values as zero correlation
  corrng[[4]][is.na(corrng[[4]])] <- 0

  # retrieve start/end of confidence intervals and median correlation coefficient
  replicates <- dim(corrng[[4]])[3]
  tail <- round(0.5*(1-confidence)*replicates)
  start_index <- tail
  end_index <- replicates - tail + 2

  res1 <- new("CorMatrix", coord1 = data.frame(corrng[[1]],corrng[[2]],corrng[[3]]), coord2 = data.frame(corrng[[1]],corrng[[2]],corrng[[3]]), cormat = corrng[[5]])
  res2 <- new("CorMatrix", coord1 = data.frame(corrng[[1]],corrng[[2]],corrng[[3]]), coord2 = data.frame(corrng[[1]],corrng[[2]],corrng[[3]]), cormat = as.data.frame(corrng[[4]][, , start_index]))
  res3 <- new("CorMatrix", coord1 = data.frame(corrng[[1]],corrng[[2]],corrng[[3]]), coord2 = data.frame(corrng[[1]],corrng[[2]],corrng[[3]]), cormat = as.data.frame(corrng[[4]][, , end_index]))
  res4 <- new("CorMatrix", coord1 = data.frame(chr,start,end), coord2 = data.frame(chr,start,end), cormat = as.data.frame(corrng[[4]][, , 1+replicates/2]))

  # return correlation coefficients and start/end positions of confidence intervals
  return(list(res1,res2,res3,res4))
}
