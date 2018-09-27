# set path for I/O
inpath <- "your_path_to_data"
promoter_path <- "your_path_to_promoter_annotations"
enhancer_path <- "your_path_to_enhancer_annotations"
outpath <- "your_path_for_matrix"

# read promoter and enhancer annotations
promoters <- readBed(promoter_path)
enhancers <- readBed(enhancer_path)

# make accessibility matrices for promoters and enhancers
amp <- makeAccMatrix(inpath, promoters)
ame <- makeAccMatrix(inpath, enhancers)

# remove empty rows
amp <- removeEmptyRows(amp)
ame <- removeEmptyRows(ame)

# get numbers of accessible promoters and enhancers
np <- dim(amp@accmat)[1]
ne <- dim(ame@accmat)[1]

# make enhancer-promoter cross correlation matrix for chromosome 1
amp1 <- cropAccMatrix(amp, 1, 1, 249250621)
ame1 <- cropAccMatrix(ame, 1, 1, 249250621)
ccor <- makeCorMatrix(amp1, ame1)

# find highest-correlated enhancer for each promoter (within a given genomic window)
window <- 100000

res <- unname(amp1[,1:3])
res[,4] <- res[,1]
res[,5] <- 0
res[,6] <- 0

for(i in 1:dim(ccor)[1]) { # loop through promoters
  # select enhancers within genomic window
  elist <- which(amp1@coord[i,1]==ame1@coord[,1] & abs((amp1@coord[i,2]+amp1@coord[i,3])/2-(ame1@coord[,2]+ame1@coord[,3])/2)<=window)

  if(length(elist)>0) {
    # find the enhancer with the highest correlation
    index <- which(ccor@cormat[i,elist] == max(ccor@cormat[i,elist]))
    if(length(index)>0) {res[i,5:6] <- ame1@coord[elist[index[1]],2:3]} # multiple enhancers found
    else {res[i,5:6] <- ame1@coord[elist[index],2:3]} # one enhancer found
  }
}

# make list of promoters with assigned enhancer
assigned <- res[res[,5]>0 & res[,6]>0,]

# save list
write.table(assigned, file=outpath, quote=F, row.names=F, col.names=F)
