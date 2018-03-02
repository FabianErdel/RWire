#' This function makes accessibility matrices from a set of bed files
#'
#' @param path Path to bed files
#' @param rois GRanges object containing ROIs
#' @param nmax Number of cells used for the analysis. When set to 0, all cells will be considered.
#' @return Accessibility matrix
#' @export

makeAccMatrix<-function(path, rois, nmax = 0) {
  # check input arguments
  if(!file.exists(path)) stop("file not found")
  if(!is.numeric(nmax)) stop("nmax must be a number")

  # determine number of ROIs
  nrois <- length(rois)
  
  # determine number of bed files in path
  files <- list.files(path, "*.bed")
  ncells <- length(files)
  
  # reduce ncells if applicable
  if((nmax > 0) & (ncells > nmax)) {ncells = nmax}

  # initialize data frame
  accmat = data.frame(chr = seqnames(rois), start = start(rois), end = end(rois))
  
  # loop through bed files
  for(i in 1:ncells) {
    # print progress
    print(paste0(i,"/",ncells,":   ",path,"/",files[i]))
    
    # read bed file
    data <- readBed(paste0(path,"/",files[i]))
    
    # count reads in ROIs
    cnt <- countReads(data, rois)
  
    # write counts in accessibility matrix
    accmat[,i+3] <- cnt
  }
  
  # add column names
  # colnames(accmat) <- c("chr", "start", "end", paste0("cell",seq(1:ncells)))
  colnames(accmat) <- c("chr", "start", "end", paste0(path,"/",files[1:ncells]))

  # return accessibility matrix
  return(accmat)
}