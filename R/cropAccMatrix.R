#' This function crops accessibility matrices
#'
#' @param am Accessibility matrix
#' @param chr Chromosome for genomic region of interest
#' @param start Start position of genomic region of interest
#' @param end End position of genomic region of interest
#' @return Accessibility matrix without empty rows
#' @export

cropAccMatrix<-function(am, chr, start, end) {
  # check input arguments
  if(class(am)!="AccMatrix") stop("am must be an AccMatrix object")
  if(!is.numeric(chr)) stop("chr must be a number")
  if(!is.numeric(start)) stop("start must be a number")
  if(!is.numeric(end)) stop("end must be a number")

  # indices for rows that should be kept
  indices <- am@coord[,1]==paste0('chr',chr) & am@coord[,2]>start & am@coord[,3]<end

  # remove empty rows
  am@accmat <- am@accmat[indices, ]
  am@coord <- am@coord[indices, ]

  # return cropped accessibility matrix
  return(am)
}
