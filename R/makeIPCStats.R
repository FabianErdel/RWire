#' Calculate statistics for the number of integrations per cell
#'
#' @param matrix Accessibility matrix
#' @export

makeIPCStats<-function(matrix) {
  # check input arguments
  if(!is.data.frame(matrix)) stop("matrix must be a data frame")
  
  # count total number of integrations for each cell
  integrations_per_cell <- colSums(matrix[,4:dim(matrix)[2]])
  
  # make and plot a histogram of integration numbers
  print(paste("Min: ",min(integrations_per_cell)))
  print(paste("Max: ",max(integrations_per_cell)))
  print(paste("Median: ",median(integrations_per_cell)))
}
