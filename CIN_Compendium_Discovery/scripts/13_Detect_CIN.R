#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

## Load args
SUMMARYFILE=args[1]
print(paste("SUMMARYFILE:", SUMMARYFILE))
DCIN=as.numeric(args[2])
print(paste("DCIN:", DCIN))
SEGMENTS=args[3]
print(paste("SEGMENTS:", SEGMENTS))
OUTPATH=args[4]
print(paste("OUTPATH:", OUTPATH))

# ## Testing
# SUMMARYFILE="/Users/drews01/data/phd/cnsigs2/data/metadata/summary.ascatTCGA.penalty70.txt"
# DCIN=19
# SEGMENTS="/Users/drews01/data/phd/cnsigs2/data/rawdata/combined.ascat.segments.smoothednormals.rds"
# OUTPATH="/Users/drews01/data/phd/cnsigs2/data/3_Pancancer_Signatures"

library(data.table)
dir.create(OUTPATH, recursive = TRUE, showWarnings = FALSE)


## Detect CIN
meta = fread(SUMMARYFILE)
meta$dCIN = meta$CNAs > DCIN
cinSamples = meta$name[ meta$dCIN ]
cinSamples = cinSamples[ ! is.na(cinSamples) ]

write.table(meta, SUMMARYFILE, sep = "\t", row.names = FALSE, quote = FALSE)

## Filter segments
segments = readRDS(SEGMENTS)
cinSegments = segments[ as.character(segments$sample) %in% cinSamples, ]
cinSegments$sample = factor(cinSegments$sample)

saveRDS(cinSegments, file.path(OUTPATH, "0_TCGA_Segments_dCIN.rds"))
write.table(cinSegments, file.path(OUTPATH, "0_TCGA_Segments_dCIN.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

