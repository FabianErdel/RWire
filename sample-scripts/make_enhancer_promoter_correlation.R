# set path for I/O
inpath <- "your_path_to_data"
promoter_path <- "your_path_to_promoter_annotations"
enhancer_path <- "your_path_to_enhancer_annotations"
outpath <- "your_path_for_matrix"

# read promoter and enhancer annotations
promoters <- readBed(promoter_path)
enhancers <- readBed(enhancer_path)

# make accessibility matrices for promoters and enhancers
amp <- makeAccMatrix2(inpath, promoters)
ame <- makeAccMatrix2(inpath, enhancers)

# remove empty rows
amp <- amp[rowSums(amp[,4:dim(amp)[2]])>0, ]
ame <- ame[rowSums(ame[,4:dim(ame)[2]])>0, ]

# binarize (optional)
amp[,4:dim(amp)[2]][amp[,4:dim(amp)[2]] > 1] <- 1
ame[,4:dim(ame)[2]][ame[,4:dim(ame)[2]] > 1] <- 1

# get numbers of accessible promoters and enhancers
np <- dim(amp)[1]
ne <- dim(ame)[1]

# make enhancer-promoter cross correlation matrix for chromosome 1
amp1 <- amp[amp[,1] == "chr1",]
ame1 <- ame[ame[,1] == "chr1",]
ccor <- makeCrossCorMatrix(amp1, ame1)

# find highest-correlated enhancer for each promoter (within a given genomic window)
window <- 100000

mat <- as.matrix(ccor[4:dim(ccor)[1],4:dim(ccor)[2]])
mat <- apply(mat, c(1,2), as.numeric)
res <- unname(amp1[,1:6])
res[,4] <- res[,1]
res[,5] <- 0
res[,6] <- 0

for(i in 1:dim(ccor)[1]) { # loop through promoters
  # select enhancers within genomic window
  elist <- which(amp1[i,1]==ame1[,1] & abs((amp1[i,2]+amp1[i,3])/2-(ame1[,2]+ame1[,3])/2)<=window)

  if(length(elist)>0) {
    # find the enhancer with the highest correlation
    index <- which(mat[i,elist] == max(mat[i,elist]))
    if(length(index)>0) {res[i,5:6] <- ame1[elist[index[1]],2:3]} # multiple enhancers found
    else {res[i,5:6] <- ame1[elist[index],2:3]} # one enhancer found
  }
}

# make list of promoters with assigned enhancer
assigned <- res[res[,5]>0 & res[,6]>0,]

# save list
write.table(assigned, file=outpath, quote=F, row.names=F, col.names=F)
