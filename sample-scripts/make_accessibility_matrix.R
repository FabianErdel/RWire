# make tiles for chromosome 18, using default settings
tiles <- makeTiles(18)

# make accessibility matrix
am <- makeAccMatrix("your_path", tiles)

# save accessibility matrix
write.table(am, file = “your_path/am", sep = "\t", row.names = FALSE, quote = FALSE)

# make histogram for integrations per cell and integrations per region
hist_ipc <- makeIPCHist(am)
print(hist_ipc)

hist_ipr <- makeIPRHist(am)
print(hist_ipr[[1]])

# print stats
makeIPCStats(am)
makeIPRStats(am)