#' This function removes rows below a certain cutoff from accessibility matrices
#'
#' @param am Accessibility matrix
#' @param cutoff Percentage of cells required to keep a row (should be between 0-100)
#' @return Accessibility matrix without rows below a certain cutoff
#' @export

removeRowsBelowCutoff<-function(am, cutoff) {
  # check input arguments
  if(class(am)!="AccMatrix") stop("am must be an AccMatrix object")
  if(!is.numeric(cutoff)) stop("cutoff must be a number")

  # indices for rows that should be kept
  indices <- rowSums(am@accmat>0)>cutoff*dim(am@accmat)[2]/100

  # remove empty rows
  am@accmat <- am@accmat[indices, ]
  am@coord <- am@coord[indices, ]

  # return accessibility matrix
  return(am)
}
