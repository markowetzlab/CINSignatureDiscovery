#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

## General configs
SOURCE=args[1]
print(paste("SOURCE:", SOURCE))
INPUTFILE=args[2]
print(paste("INPUTFILE:", INPUTFILE))
OUTPATH=args[3]
print(paste("OUTPATH:", OUTPATH))
CORES=as.numeric(args[4])
print(paste("CORES:", CORES))

# Step 1: Extract copy-number features
PREPATH=args[5]
print(paste("PREPATH:", PREPATH))
RMNORM=as.logical(args[6])
print(paste("RMNORM:", RMNORM))

### Testing
## General configs
# BASE="/Users/drews01/data/phd/cnsigs2/data"
# SOURCE="/Users/drews01/data/phd/cnsigs2/finalScripts/00_main_functions.R"
# INPUTFILE=file.path(BASE, "3_TCGA_specific_signatures/0_allSegments_filtered_purity0.4.rds")
# OUTPATH=file.path(BASE, "3_TCGA_specific_signatures")
# CORES=7
# 
# # Step 1: Extract copy-number features
# PREPATH="/Users/drews01/data/phd/cnsigs2/data/Geoffs_data/"
# RMNORM=TRUE

## 0) Prepare run
dir.create(OUTPATH, recursive = TRUE)
source(SOURCE)

## 1) extract copy-number features from unrounded ASCAT calls
out1 = file.path(OUTPATH, "1_tcga_filtered_ecnf.rds")
if( ! file.exists(out1) ) {
    tcga.table = readRDS(INPUTFILE)
    # Set factors again so no empty samples are present. To be sure.
    tcga.table$sample = factor(tcga.table$sample)
    tcga.list = split( tcga.table, tcga.table$sample )
    tcga.ecnf = extractCopynumberFeatures( tcga.list, cores = CORES, prePath=PREPATH, rmNorm = RMNORM )
    saveRDS(tcga.ecnf, file = out1 )
} else {
    tcga.ecnf = readRDS( out1 )
}

