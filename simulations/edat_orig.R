library(curatedMetagenomicData)

studies = unique(combined_metadata$dataset_name[which((!combined_metadata$BMI %in% NA | !combined_metadata$cholesterol %in% NA) & combined_metadata$body_site %in% "stool")])

metadata = curatedMetagenomicData(paste0(studies, ".marker_abundance.stool"), dryrun = F) 

eset_orig = vector("list", length(metadata))
for (i in 1:length(metadata)) {
  eset_orig[[i]] = t(exprs(metadata[[i]]))[, 1:10000]
}

# Work with the intersection of rows
cn <- lapply(eset_orig, colnames)
cn_int <- Reduce(intersect, cn)

edat <- vector("list", length(metadata))

for(i in 1:length(edat)){
  edat[[i]] <- eset_orig[[i]][,cn_int]
}

edat_orig <- edat

# Normalize the columns
for(i in 1:length(edat_orig)){
  edat_orig[[i]] <- apply(edat_orig[[i]], 2, scale)
}


save(edat_orig, file="edat_orig.RData")


