#' This function saves a matrix as an image
#'
#' @param mat Matrix
#' @param path Storage path for image
#' @param w Image width
#' @param h Image height
#' @param color Color LUT to be used
#' @param min Minimum value for color LUT
#' @param max Maximum value for color LUT
#' @export

saveImage<-function(mat, path, w, h, color, min = 0, max = -1) {
  # check input arguments
  if(!is.matrix(mat)) stop("mat must be a matrix")
  if(!is.character(path)) stop("path must be a character string")

  # assign max value
  if(max < min) {max = max(mat)}
  
  # save matrix
  png(path, width = w, height = h)
  par(mar = rep(0, 4))
  image(mat, zlim=c(min, max), axes = FALSE, col = color)
  dev.off()
}