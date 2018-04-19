# set path for I/O
inpath <- "your_path_to_data"
outpath1 <- "your_path_for_matrix"
outpath2 <- "your_path_for_matrix_image"
outpath3 <- "your_path_for_scaled_matrix"

# read accessibility matrix
am <- fread(inpath, data.table=FALSE, sep="\t")

# remove empty rows
am <- am[rowSums(am[,4:dim(am)[2]])>0, ]

# select genomic region of interest
start <- 1
end <- 1000000

# crop accessibility matrix
am <- am[am[,2]>start & am[,3]<end, ]

# get number of accessible regions
nrois <- dim(am)[1]-3

# make correlation matrix
cor <- makeCorMatrix(am)

# save correlation matrix (as table)
write.table(cor, file = outpath1, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# save correlation matrix (as image)
saveImage(as.matrix(cor[,4:(nrois+3)]), outpath2, nrois, nrois, grey(seq(0, 1, length = 256)), -1, 1)

# make correlation matrix scaled according to genomic positions
# each pixel will corresponds to a genomic range of (end-start)/size
size <- 200

# make scaled matrix
chr <- 1
scm <- scaleCorMatrix(cor, chr, size)

# save matrix
saveImage(scm, outpath3, size, size, getCorLUT(), -1, 1)