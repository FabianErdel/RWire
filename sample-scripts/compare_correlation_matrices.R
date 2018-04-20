# set path for I/O
inpath1 <- "your_path_to_data"
inpath2 <- "your_path_to_data"
outpath <- "your_path_for_matrix"

# read accessibility matrices (containing data for the same genomic regions, but for different conditions)
ram1 <- fread(inpath1, data.table=FALSE, sep="\t")
ram2 <- fread(inpath2, data.table=FALSE, sep="\t")

# remove empty rows
am1 <- ram1[rowSums(ram1[,4:dim(ram1)[2]])>0 | rowSums(ram2[,4:dim(ram2)[2]])>0, ]
am2 <- ram2[rowSums(ram1[,4:dim(ram1)[2]])>0 | rowSums(ram2[,4:dim(ram2)[2]])>0, ]

# get number of accessible regions
nrois <- dim(am1)[1]

# make correlation matrices
cor1 <- makeCorMatrix(am1)
cor2 <- makeCorMatrix(am2)

# get bootstrapped correlation coefficients
replicates <- 100
cr1 <- getCorRanges(am1, replicates)
cr2 <- getCorRanges(am2, replicates)

# compare bootstrapped correlation coefficients and calculate p-values to assess their difference (KS test)
pval <- compareCorrelations(cr1, cr2)

# set p-values for correlations that are zero in both matrices to unity
pval[cor1==0 & cor2==0] <- 1

# save matrix with p-values as image
saveImage(pval, outpath, nrois, nrois, getCorLUT(), -1, 1)
