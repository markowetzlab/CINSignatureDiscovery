#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

BASE = args[1]
print(paste("BASE:", BASE))
INFOLDER = args[2]
print(paste("INFOLDER:", INFOLDER))
OUTFOLDER = args[3]
print(paste("OUTFOLDER:", OUTFOLDER))
CHRSIZES = args[4]
print(paste("CHRSIZES:", CHRSIZES))
FUNSCRIPT = args[5]
print(paste("FUNSCRIPT:", FUNSCRIPT))
CORES = as.numeric(args[6])
print(paste("CORES:", CORES))

## Testing
# BASE = "/Users/drews01/data/phd/cnsigs2/data"
# OUTFOLDER = "0_cnprofiles/ascat"
# INFOLDER = "rawdata/split-by-cancer-type"
# CHRSIZES = file.path(BASE, "Geoffs_data/hg19.chrom.sizes.txt")
# FUNSCRIPT = file.path(BASE, "../scripts/0_cnprofiles_functions.R")
# CORES=10

dir.create(file.path(BASE, OUTFOLDER), showWarnings = FALSE, recursive = TRUE)
setwd(BASE)
source(FUNSCRIPT)
library(foreach)
library(doMC)

# plot per sample and per cancer types for all listed files
# Simple:
# plotMaster(BASE, INFOLDER, OUTFOLDER, CHRSIZES)

# Advanced:
chrSizes = getChrSizes(CHRSIZES)
cancerFiles = list.files(file.path(BASE, INFOLDER), pattern = "*.filt.*")
registerDoMC(cores=CORES)
foreach(thisCancerFile=cancerFiles) %dopar% {
    
    thisCancer = unlist(strsplit(thisCancerFile, "[.]"))[1]
    print(thisCancer)
    dir.create(file.path(BASE, OUTFOLDER, thisCancer), recursive = TRUE, showWarnings = FALSE)
    
    thisFileName = file.path(BASE, INFOLDER, thisCancerFile)
    eventsPerSample = loadEvents( thisFileName )
    
    # Plot per sample
    pSamples = plotPerSamples(eventsPerSample, BASE, OUTFOLDER, thisCancer, chrSizes)
    
    # Cancer overview pdf
    plotsPerCancer(pSamples, BASE, OUTFOLDER, thisCancer)
    
}

## With rounded values, ABSOLUTE data etc. see script 0_cnprofiles.R in scripts
