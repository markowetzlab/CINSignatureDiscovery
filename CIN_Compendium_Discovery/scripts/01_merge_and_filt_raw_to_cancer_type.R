#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

library(foreach)

## Testing
# THISPATH = "~/phd/prjcts/cnsigs2/data/rawdata"
# INPUT = file.path( THISPATH, "rawsegments" )
# OUTPUTSEGMENTS = file.path( THISPATH, "split-by-cancer-type" )
# OUTPUTMETA = file.path( THISPATH, "../metadata" )
# dir.create(OUTPUT, recursive = TRUE)
# SUMMARYFILE = file.path( THISPATH, "summary.ascatTCGA.penalty70.txt" )
# PURITY = 0.4

## Transer variables from command line options
THISPATH = args[1]
print(paste("THISPATH:", THISPATH))
INPUT = file.path( THISPATH, args[2] )
print(paste("INPUT:", INPUT))
SUMMARYFILE = file.path( THISPATH, args[3] )
print(paste("SUMMARYFILE:", SUMMARYFILE))

OUTPUTSEGMENTS = file.path( THISPATH, args[4] )
print(paste("OUTPUTSEGMENTS:", OUTPUTSEGMENTS))
OUTPUTMETA = file.path( THISPATH, args[5] )
print(paste("OUTPUTMETA:", OUTPUTMETA))

PURITY = as.numeric(args[6])
print(paste("PURITY:", PURITY))


## Load
dir.create(OUTPUTSEGMENTS, recursive = TRUE)
dir.create(OUTPUTMETA, recursive = TRUE)

dfSummary = read.table( SUMMARYFILE, header = TRUE )
allFiles = list.files(INPUT, pattern = "TCGA")

# samples have to have an ASCAT solution ("solution"), shall be the representative for their sample in case of multisamples ("rep"),
# pass Kerstin's QC ("pass"), have a purity less than 100% (100% is most likely an ASCAT artefact) and not be oversegmented due to
# a mismatch between logR and BAF values (some samples have wild discordances there).
dfSingleSamples = dfSummary[ dfSummary$rep == TRUE & dfSummary$solution == TRUE & dfSummary$pass == TRUE &
                                 dfSummary$purity < 1 & dfSummary$notOverSegm == TRUE & dfSummary$purity > PURITY, ]

cancerTypes = unique( dfSingleSamples$cancer_type )
for(thisCancer in cancerTypes) {

    print(thisCancer)
    theseSamples = dfSingleSamples[ dfSingleSamples$cancer_type == thisCancer, ]$name
    theseFiles = allFiles[ allFiles %in% paste0(theseSamples, ".segments.raw.txt") ]
    theseSegments = foreach( thisFile=theseFiles, .combine='rbind' ) %do% {
        return( read.table( file.path(INPUT, thisFile), header = TRUE) )
    }
    write.table(theseSegments, file = file.path(OUTPUTSEGMENTS, paste0(thisCancer, ".ascat.segments.filt.txt")),
                quote=FALSE, row.names=FALSE, sep="\t")

}
