# set path for I/O
inpath1 <- "your_path_to_data1"
inpath2 <- "your_path_to_data2"
outpath1 <- "your_path_for_matrix"
outpath2 <- "your_path_for_matrix_image"

# read accessibility matrices
am1 <- data.table::fread(inpath1, data.table=FALSE, sep="\t")
am2 <- data.table::fread(inpath2, data.table=FALSE, sep="\t")

# remove empty rows
am1 <- removeEmptyRows(am1)
am2 <- removeEmptyRows(am2)

# get numbers of accessible regions
nrois1 <- dim(am1@accmat)[1]
nrois2 <- dim(am2@accmat)[1]

# make cross correlation matrix
ccor <- makeCorMatrix(am1, am2)

# save cross correlation matrix (as table)
write.table(ccor@cormat, file = outpath1, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# save correlation matrix (as image)
saveImage(as.matrix(ccor@cormat), outpath2, nrois1, nrois2, grey(seq(0, 1, length = 256)), -1, 1)
