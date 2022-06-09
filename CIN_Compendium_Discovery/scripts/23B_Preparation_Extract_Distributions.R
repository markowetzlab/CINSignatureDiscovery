#!/usr/bin/env Rscript

## Extract segment size and changepoint distributions from TCGA ECNF to be modelled with VI DPGMM 

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

# General
INPUT=args[1]
print(paste("INPUT:", INPUT))
OUTPUT=args[2]
dir.create(OUTPUT, showWarnings = FALSE, recursive = TRUE)
print(paste("OUTPUT:", OUTPUT))

# ## Testing
# INPUT="~/SignatureDiscovery/CIN_Compendium_Discovery/data/backup/1_tcga_filtered_ecnf.rds"
# OUTPUT="~/SignatureDiscovery/CIN_Compendium_Discovery/data/backup"

# Read in
lIn = readRDS(INPUT)

# Extract segment size
dfSegsize = lIn$segsize
dfSegsize$value

write.table(dfSegsize$value, file.path(OUTPUT, "1_tcga_segmentsize_dist.csv"), quote = FALSE, sep=",", row.names = FALSE, col.names = FALSE)

# Extract changepoint
dfChangepoint = lIn$changepoint
dfChangepoint$value

write.table(dfChangepoint$value, file.path(OUTPUT, "1_tcga_changepoint_dist.csv"), quote = FALSE, sep=",", row.names = FALSE, col.names = FALSE)
