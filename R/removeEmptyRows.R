#' This function removes empty rows from accessibility matrices
#'
#' @param am Accessibility matrix
#' @return Accessibility matrix without empty rows
#' @export

removeEmptyRows<-function(am) {
  # check input arguments
  if(class(am)!="AccMatrix") stop("am must be an AccMatrix object")

  # indices for rows that should be kept
  indices <- rowSums(am@accmat)>0

  # remove empty rows
  am@accmat <- am@accmat[indices, ]
  am@coord <- am@coord[indices, ]

  # return accessibility matrix
  return(am)
}
