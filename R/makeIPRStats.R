#' Calculate statistics for the number of integrations per ROI
#'
#' @param am Accessibility matrix
#' @export

makeIPRStats<-function(am) {
  # check input arguments
  if(class(am)!="AccMatrix") stop("am must be an AccMatrix object")

  # count total number of integrations for each ROI
  integrations_per_roi <- rowSums(am@accmat)

  # make and plot a histogram of integration numbers
  print("Integration statistics")
  print(paste("Min: ",min(integrations_per_roi)))
  print(paste("Max: ",max(integrations_per_roi)))
  print(paste("Median: ",median(integrations_per_roi)))

  # count total number of cells with integrations for each ROI
  cells_per_roi <- rowSums(am@accmat > 0)

  # make and plot a histogram of integration numbers
  print("Cells with integration statistics")
  print(paste("Min: ",min(cells_per_roi)))
  print(paste("Max: ",max(cells_per_roi)))
  print(paste("Median: ",median(cells_per_roi)))
}
