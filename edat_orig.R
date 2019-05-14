
# the code below is from https://github.com/prpatil/csml_rep/blob/master/csml_runner.R

#data(package="curatedOvarianData")
source(system.file("extdata", "patientselection.config",package="curatedOvarianData"))
sapply(ls(), function(x) if(!x %in% c("remove.samples", "duplicates")) print(get(x)))
source(system.file("extdata", "createEsetList.R", package = "curatedOvarianData"))
# Remove esets with missing gene expression data
ridx <- which(unlist(lapply(esets, function(x){sum(is.na(exprs(x)))})) > 0)
esets <- esets[-ridx]
eset_orig <- vector("list", length(esets))

# Convert eset list to set of matrices
for(i in 1:length(esets)){
  eset_orig[[i]] <- t(exprs(esets[[i]]))
}

# Work with the intersection of rows
cn <- lapply(eset_orig, colnames)
cn_int <- Reduce(intersect, cn)

edat <- vector("list", length(esets))

for(i in 1:length(edat)){
  edat[[i]] <- eset_orig[[i]][,cn_int]
}

edat_orig <- edat

# Normalize the columns
for(i in 1:length(edat_orig)){
  edat_orig[[i]] <- apply(edat_orig[[i]], 2, scale)
}


save(edat_orig, file="edat_orig.RData")


