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
w <- 100000 # 100 kb window
pepairs <- getPEPairs(ccor, window = w, max.only = T)

# make list of promoters with assigned enhancer
assigned <- pepairs@promoters[,1:3]
assigned[,4:6] <- pepairs@enhancers[pepairs@promoters$`assigned enhancer indices`,1:3]

# save list
write.table(assigned, file=outpath, quote=F, row.names=F, col.names=F)
