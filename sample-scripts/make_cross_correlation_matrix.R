# set path for I/O
inpath1 <- "your_path_to_data1"
inpath2 <- "your_path_to_data2"
outpath1 <- "your_path_for_matrix"
outpath2 <- "your_path_for_matrix_image"

# read accessibility matrices
am1 <- fread(inpath1, data.table=FALSE, sep="\t")
am2 <- fread(inpath2, data.table=FALSE, sep="\t")

# remove empty rows
am1 <- am1[rowSums(am1[,4:dim(am1)[2]])>0, ]
am2 <- am2[rowSums(am2[,4:dim(am2)[2]])>0, ]

# get numbers of accessible regions
nrois1 <- dim(am1)[1]
nrois2 <- dim(am2)[1]

# make cross correlation matrix
ccor <- makeCrossCorMatrix(am1, am2)

# save cross correlation matrix (as table)
write.table(ccor, file = outpath1, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# save correlation matrix (as image)
saveImage(data.matrix(ccor[4:(nrois1+3), 4:(nrois2+3)]), outpath2, nrois1, nrois2, grey(seq(0, 1, length = 256)), -1, 1)
