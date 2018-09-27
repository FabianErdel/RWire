# set path for I/O
inpath1 <- "your_path_to_data1"
inpath2 <- "your_path_to_data2"
outpath1 <- "your_path_for_matrix"
outpath2 <- "your_path_for_scaled_matrix"

# read accessibility matrices (containing data for the same genomic regions, but for different conditions)
ram1 <- data.table::fread(inpath1, data.table=FALSE, sep="\t")
ram2 <- data.table::fread(inpath2, data.table=FALSE, sep="\t")

# remove empty rows
am1 <- removeEmptyRows(ram1)
am2 <- removeEmptyRows(ram2)

# get number of accessible regions
nrois <- dim(am1@accmat)[1]

# make correlation matrices
cor1 <- makeCorMatrix(am1)
cor2 <- makeCorMatrix(am2)

# get bootstrapped correlation coefficients
replicates <- 100
cr1 <- getCorRanges(am1, replicates)
cr2 <- getCorRanges(am2, replicates)

# calculate p-values for the difference between correlation coefficients (based on bootstrap distributions)
pval <- compareCorCoeffs(cr1, cr2)

# set p-values for correlations that are zero in both matrices to unity
pval[[4]][cor1@cormat==0 & cor2@cormat==0] <- 1

# make matrix with exponents of p-values
logpval <- -log10(pval[[4]])
maxval <- max(logpval[is.finite(logpval)])
logpval[!is.finite(logpval)] <- maxval

# save p-value matrix as image
saveImage(logpval, outpath1, nrois, nrois, getCorLUT(), min(logpval), max(logpval))

# make p-value matrix scaled according to genomic positions (for chromosome 1)
size <- 200
spval <- scaleCorMatrix(cbind(pval[[1]], pval[[2]], pval[[3]], data.frame(logpval)), chr=1, size=size)

# save scaled p-value matrix as image
saveImage(spval, outpath2, size, size, getCorLUT(), min(spval), max(spval))
