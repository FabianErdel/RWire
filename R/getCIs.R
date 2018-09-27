#' This function returns confidence intervals
#'
#' @param corrng Table with annotations and correlation coefficients
#' @param confidence Confidence level (between 0 and 1)
#' @return Stack of CorMatrices containing correlation coefficients (first layer) along with minimum, maximum and median values from bootstrapping (second, third and fourth layer)
#' @export

getCIs<-function(corrng, confidence = 0.95) {
  # check input arguments
  if(!is.list(corrng)) stop("corrng must be a list")

  # retrieve items from corrng list
  chr <- corrng[[1]]
  start <- corrng[[2]]
  end <- corrng[[3]]
  cortbl <- corrng[[4]]


  # retrieve start/end of confidence intervals and median correlation coefficient
  replicates <- dim(cortbl)[3]
  tail <- round(0.5*(1-confidence)*replicates)
  start_index <- tail
  end_index <- replicates - tail + 2

  res1 <- new("CorMatrix", coord1 = cbind(chr,start,end), coord2 = cbind(chr,start,end), cormat = corrng[[5]])
  res2 <- new("CorMatrix", coord1 = cbind(chr,start,end), coord2 = cbind(chr,start,end), cormat = as.data.frame(cortbl[, , start_index]))
  res3 <- new("CorMatrix", coord1 = cbind(chr,start,end), coord2 = cbind(chr,start,end), cormat = as.data.frame(cortbl[, , end_index]))
  res4 <- new("CorMatrix", coord1 = cbind(chr,start,end), coord2 = cbind(chr,start,end), cormat = as.data.frame(cortbl[, , 1+replicates/2]))

  # return correlation coefficients and start/end positions of confidence intervals
  return(list(res1,res2,res3,res4))
}
