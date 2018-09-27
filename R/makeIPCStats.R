#' Calculate statistics for the number of integrations per cell
#'
#' @param am Accessibility matrix
#' @export

makeIPCStats<-function(am) {
  # check input arguments
  if(class(am)!="AccMatrix") stop("am must be an AccMatrix object")

  # count total number of integrations for each cell
  integrations_per_cell <- colSums(am@accmat)

  # make and plot a histogram of integration numbers
  print(paste("Min: ",min(integrations_per_cell)))
  print(paste("Max: ",max(integrations_per_cell)))
  print(paste("Median: ",median(integrations_per_cell)))
}
