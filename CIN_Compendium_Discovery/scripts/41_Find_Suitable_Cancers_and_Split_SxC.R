# Determine which cohorts and samples are suitable for a cancer-specific signature derivation

library(data.table)

args=commandArgs(trailingOnly = TRUE)

METADATA=args[1]
print(paste("METADATA:",METADATA))
MINSAMPS=as.numeric(args[2])
print(paste("MINSAMPS:",MINSAMPS))
IGNOREBRCASPLIT=as.logical(args[3])
print(paste("IGNOREBRCASPLIT:",IGNOREBRCASPLIT))
SXCMATRIX=args[4]
print(paste("SXCMATRIX:",SXCMATRIX))
OUTPUTDIR=args[5]
print(paste("OUTPUTDIR:",OUTPUTDIR))

# ## Testing
# METADATA="/Users/drews01/phd/prjcts/cnsigs2/data/metadata/summary.ascatTCGA.penalty70.txt"
# MINSAMPS=100
# IGNOREBRCASPLIT=TRUE
# SXCMATRIX="/Users/drews01/phd/prjcts/cnsigs2/data/3_Pancancer_Signatures/3_SxC_uninfPrior.rds"
# OUTPUTDIR="/Users/drews01/phd/prjcts/cnsigs2/data/4_Cancer_specific_Signatures"

metaRaw = fread(METADATA)
meta = metaRaw[metaRaw$dCIN, ]
matSxC = readRDS(SXCMATRIX)

if(IGNOREBRCASPLIT) {
    meta$cancer_type[ grepl("BRCA", meta$cancer_type) ] = "BRCA"
}

samplePerCancer = sort(table(meta$cancer_type), decreasing = TRUE)
suitableCancers = names( samplePerCancer[ samplePerCancer >= MINSAMPS ] )

# Loop over cancer types, extract matrix and save matrix as CxS output
lOut = list()
for(thisCancer in suitableCancers) {
    
    thisSxC = matSxC[ rownames(matSxC) %in% meta$name[ meta$cancer_type == thisCancer ], ]
    dir.create(file.path(OUTPUTDIR, thisCancer), showWarnings = FALSE, recursive = TRUE)
    saveRDS(thisSxC, file.path(OUTPUTDIR, thisCancer, paste0("3_SxC_", thisCancer, ".rds")))
    write.table(thisSxC, file = file.path(OUTPUTDIR, thisCancer, paste0("3_SxC_", thisCancer, ".txt")), quote = FALSE, sep = "\t", 
                row.names = TRUE, col.names = TRUE)
    thisCxS = t( thisSxC )
    saveRDS(thisCxS, file.path(OUTPUTDIR, thisCancer, paste0("3_CxS_", thisCancer, ".rds")))
    write.table(thisCxS, file = file.path(OUTPUTDIR, thisCancer, paste0("3_CxS_", thisCancer, ".txt")), quote = FALSE, sep = "\t", 
                row.names = TRUE, col.names = TRUE)
    
    # Just save number of patients and cancer for methods section
    lOut[[length(lOut)+1]] = c(thisCancer, dim(thisSxC)[1])
    
}

dtOut = data.table(do.call(rbind, lOut))
dtOut = dtOut[ order(dtOut$V2, decreasing = TRUE), ]
colnames(dtOut) = c("Cancer", "SampleNumber")
write.table(dtOut, file = file.path(OUTPUTDIR, paste0("Overview_Cancer_specific_data_sets.txt")), quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)
