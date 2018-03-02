#' This function makes a histogram for the number of integrations per ROI
#'
#' @param matrix Accessibility matrix
#' @return Histogram (as a plot object)
#' @export

makeIPRHist<-function(matrix) {
  # check input arguments
  if(!is.data.frame(matrix)) stop("matrix must be a data frame")
  
  # count total number of cells with integrations for each roi
  cells_per_roi <- rowSums(matrix[,4:dim(matrix)[2]] > 0) # count columns with values above zero
  
  # make, fit and plot a histogram of cell numbers with integrations
  
  # fit negative binomial distribution
  clspr_fit_nbinom <- fitdistrplus::fitdist(cells_per_roi[(cells_per_roi > 0) & (cells_per_roi < dim(matrix)[2]/10)], "nbinom")
  
  cutoff <- clspr_fit_nbinom$estimate[[2]] + 4*sqrt(clspr_fit_nbinom$estimate[[2]]*(1+clspr_fit_nbinom$estimate[[2]]/clspr_fit_nbinom$estimate[[1]]))
  
  plot <- ggplot2::qplot(cells_per_roi, binwidth = 1) + 
          geom_histogram(binwidth = 1) + scale_y_log10(limits=c(1, NA)) + stat_function(fun = function(x, size, mu, n) { 
          ((ceiling(x)-x)*dnbinom(x = floor(x), size = size, mu = mu) + (x-floor(x))*dnbinom(x = ceiling(x), size = size, mu = mu) + (1-ceiling(x)+floor(x))*dnbinom(x = floor(x), size = size, mu = mu)) * n}, args = c(size = as.numeric(clspr_fit_nbinom$estimate[1]), mu = as.numeric(clspr_fit_nbinom$estimate[2]), n = dim(matrix)[1]), col = 'blue') +
          labs(title = "Number of 'cells with integrations' per region", x = "Cells with integrations", y = "Regions") +
          geom_vline(xintercept = cutoff, linetype="dotted")

  return(list(plot, cutoff))
}