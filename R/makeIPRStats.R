#' Calculate statistics for the number of integrations per ROI
#'
#' @param matrix Accessibility matrix
#' @export

makeIPRStats<-function(matrix) {
  # check input arguments
  if(!is.data.frame(matrix)) stop("matrix must be a data frame")
  
  # count total number of integrations for each ROI
  integrations_per_roi <- rowSums(matrix[,4:dim(matrix)[2]])
  
  # make and plot a histogram of integration numbers
  print("Integration statistics")
  print(paste("Min: ",min(integrations_per_roi)))
  print(paste("Max: ",max(integrations_per_roi)))
  print(paste("Median: ",median(integrations_per_roi)))
  
  # count total number of cells with integrations for each ROI
  cells_per_roi <- rowSums(matrix[,4:dim(matrix)[2]] > 0)
  
  # make and plot a histogram of integration numbers
  print("Cells with integration statistics")
  print(paste("Min: ",min(cells_per_roi)))
  print(paste("Max: ",max(cells_per_roi)))
  print(paste("Median: ",median(cells_per_roi)))
}