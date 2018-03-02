#' This function makes a color LUT for correlation matrices
#'
#' @param steps Number of steps for the color LUT
#' @return Color LUT
#' @export

getCorLUT<-function(steps = 255) {
  # check input arguments
  if(!is.numeric(steps)) stop("steps must be a number")

  # make vectors with steps
  steps1 <- seq(0, 127, length = steps/2)
  steps2 <- seq(128, 255, length = steps/2)
  
  # blue-white gradient for negative values
  bw <- rgb(1-(128-steps1)/128, 1-(128-steps1)/128, 1)
  
  # red-white gradient for positive values
  rw <- rgb(1, 1-(steps2-128)/128, 1-(steps2-128)/128)
  
  # return color LUT
  return(c(bw, rw))
}