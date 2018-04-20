# set path for I/O
inpath1 <- "your_path_to_data"
inpath2 <- "your_path_to_data"
outpath <- "your_path_for_matrix"


# read accessibility matrices (containing data for the same genomic regions, but for different conditions)
am1 <- fread(inpath1, data.table=FALSE, sep="\t")
am2 <- fread(inpath2, data.table=FALSE, sep="\t")

# remove empty rows
am1 <- am1[rowSums(am1[,4:dim(am1)[2]])>0, ]
am2 <- am2[rowSums(am2[,4:dim(am2)[2]])>0, ]

# select genomic region of interest
start <- 1
end <- 1000000

# crop accessibility matrices
am1 <- am1[am1[,2]>start & am1[,3]<end, ]
am2 <- am2[am2[,2]>start & am2[,3]<end, ]

# get number of accessible regions
nrois <- dim(am1)[1]-3

# get bootstrapped correlation coefficients
replicates <- 100
cr1 <- getCorRanges(am1, replicates)
cr2 <- getCorRanges(am2, replicates)

# compare bootstrapped correlation coefficients and calculate p-values to assess their difference (KS test)
pval <- compareCorrelations(cr1, cr2)

# save matrix with p-values as image
saveImage(pval, outpath, nrois, nrois, getCorLUT(), -1, 1)
