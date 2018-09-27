#' This function makes a histogram for the number of integrations per cell
#'
#' @param am Accessibility matrix
#' @return Histogram (plot object)
#' @export

makeIPCHist<-function(am) {
  # check input arguments
  if(class(am)!="AccMatrix") stop("am must be an AccMatrix object")

  # count total number of integrations for each cell
  integrations_per_cell <- colSums(am@accmat)

  # make and plot a histogram of integration numbers
  intpc_hist <- hist(integrations_per_cell, breaks="FD", plot = FALSE)
  plot <- ggplot2::qplot(integrations_per_cell, breaks = intpc_hist$breaks) + geom_histogram() +
    labs(title = "Histogram of integrations per cell", x = "Integrations", y = "Cells")

  return(plot)
}
