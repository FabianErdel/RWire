#' This function crops correlation matrices
#'
#' @param cm Correlation matrix
#' @param chr Chromosome for genomic region of interest
#' @param start Start position of genomic region of interest
#' @param end End position of genomic region of interest
#' @return Accessibility matrix without empty rows
#' @export

cropCorMatrix<-function(cm, chr, start, end) {
  # check input arguments
  if(class(cm)!="CorMatrix") stop("cm must be a CorMatrix object")
  if(!is.numeric(chr)) stop("chr must be a number")
  if(!is.numeric(start)) stop("start must be a number")
  if(!is.numeric(end)) stop("end must be a number")

  # indices for rows/columns that should be kept
  indices1 <- cm@coord1[,1]==paste0('chr',chr) & cm@coord1[,2]>start & cm@coord1[,3]<end
  indices2 <- cm@coord2[,1]==paste0('chr',chr) & cm@coord2[,2]>start & cm@coord2[,3]<end

  # remove empty rows
  cm@coord1 <- cm@coord1[indices1, ]
  cm@coord2 <- cm@coord2[indices2, ]
  cm@cormat <- cm@cormat[indices1, indices2]

  # return cropped correlation matrix
  return(cm)
}
