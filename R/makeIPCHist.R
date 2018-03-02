#' This function makes a histogram for the number of integrations per cell
#'
#' @param matrix Accessibility matrix
#' @return Histogram (as a plot object)
#' @export

makeIPCHist<-function(matrix) {
  # check input arguments
  if(!is.data.frame(matrix)) stop("matrix must be a data frame")

  # count total number of integrations for each cell
  integrations_per_cell <- colSums(matrix[,4:dim(matrix)[2]])

  # make and plot a histogram of integration numbers
  intpc_hist <- hist(integrations_per_cell, breaks="FD", plot = FALSE)
  plot <- ggplot2::qplot(integrations_per_cell, breaks = intpc_hist$breaks) + geom_histogram()

  return(plot)
}
