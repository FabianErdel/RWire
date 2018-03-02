#' This function makes a color LUT for accessibility matrices
#'
#' @param steps Number of steps for the color LUT
#' @return Color LUT
#' @export

getAccLUT<-function(steps = 255) {
  # check input arguments
  if(!is.numeric(steps)) stop("steps must be a number")

  # make vector with steps
  steps <- seq(0, 255, length = steps)

  # red-white gradient for positive values
  rw <- rgb(1, 1-steps/255, 1-steps/255)
  
  # return color LUT
  return(rw)
}