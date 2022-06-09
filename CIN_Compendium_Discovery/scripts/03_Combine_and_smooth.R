#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

## Load args
INFOLDER=args[1]
print(paste("INFOLDER:", INFOLDER))
OUTCOMBFILE=args[2]
print(paste("OUTCOMBFILE:", OUTCOMBFILE))
OUTSMOOTHFILE=args[3]
print(paste("OUTSMOOTHFILE:", OUTSMOOTHFILE))
SUMMARYFILE=args[4]
print(paste("SUMMARYFILE:", SUMMARYFILE))
WIGGLE=as.numeric(args[5])
print(paste("WIGGLE:", WIGGLE))
IGNOREDELS=as.logical(args[6])
print(paste("IGNOREDELS:", IGNOREDELS))
CORES=as.numeric(args[7])
print(paste("CORES:", CORES))

### Testing purposes
# DATA="/Users/drews01/phd/prjcts/cnsigs2/data"
# INFOLDER=file.path(DATA, "rawdata/split-by-cancer-type")
# OUTCOMBFILE=file.path(DATA, "rawdata/combined.ascat.segments.filt")
# OUTSMOOTHFILE=file.path(DATA, "rawdata/combined.ascat.segments.smoothednormals")
# SUMMARYFILE=file.path(DATA, "metadata/summary.ascatTCGA.penalty70.txt")
# WIGGLE=0.1
# IGNOREDELS=FALSE
# CORES=7

library(data.table)
library(GenomicRanges)
library(foreach)
library(doMC)

### PART 1: Combine files
## Load input
cancerFiles = list.files(INFOLDER, pattern = "*.filt.*", full.names = TRUE)
allSegments = data.frame()
lAll = lapply(cancerFiles, function(thisCancerFile) {
   
    print(thisCancerFile)
    thisCancer = read.table(thisCancerFile, header = TRUE )
    outCancer = data.table( "chromosome" = thisCancer[,2], 
                          "start" = as.numeric( thisCancer[,3] ), 
                          "end" = as.numeric( thisCancer[,4] ), 
                          "segVal" = as.numeric( thisCancer[,7] ) + as.numeric( thisCancer[,8] ), 
                          "sample" = thisCancer[,1] )
    return( outCancer )
    
} )

dtDat = rbindlist(lAll)

# explicit conversion to numeric. Just to be on the safe site
dtDat$start = as.numeric( dtDat$start )
dtDat$end = as.numeric( dtDat$end )
dtDat$segVal = as.numeric( dtDat$segVal )

## Write txt and rds
write.table(dtDat, paste0(OUTCOMBFILE, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
saveRDS(dtDat, paste0(OUTCOMBFILE, ".rds") )


### PART II: Smooth
idSmoothingTargets = function(dtDat, WIGGLE, colNameSegVal, colNameChr, IGNOREDELS = TRUE) {
    
    ### Check column name
    testSegVal = dtDat[[colNameSegVal]][1]
    testChr = dtDat[[colNameChr]][1]
    if(! is.numeric(testSegVal)) { stop("Segment Value column has no numeric value in it. Supplied correct column name?")}
    if(is.null(testSegVal)) { stop("Chromosome column has no numeric value in it. Supplied correct column name?")}
    
    # take differences to segment down below
    dtDat$diffs = c( abs( dtDat[[colNameSegVal]][1:(nrow(dtDat)-1)] - dtDat[[colNameSegVal]][2:nrow(dtDat)] ), WIGGLE+1)
    # set TRUE if difference to next segment is smaller than the user supplied cutoff
    dtDat$smooth = dtDat$diffs <= WIGGLE
    # set all segments which are last in a chromosome to FALSE. This also prevents leaking to other samples and cohorts.
    dtDat$smooth[ cumsum( rle(as.character(dtDat[[colNameChr]]))$lengths ) ] = FALSE
    
    # Ignore deletions if wished
    if(IGNOREDELS) { dtDat$smooth[ dtDat[[colNameSegVal]] == 0 ] = FALSE }
    
    return( dtDat )    
}


smoothSegments = function(lRaw, CORES, SMOOTHINGFACTOR, colNameMerge, colNameChr, colNameStart, colNameEnd, 
                          IGNOREDELS = TRUE, asDf = FALSE) {
    
    ### Check column names
    test = lRaw[[1]]
    testMerge = test[[colNameMerge]][1]
    testChr = test[[colNameChr]][1]
    testStart = test[[colNameStart]][1]
    testEnd = test[[colNameEnd]][1]
    if(! is.numeric(testMerge)) { stop("Merge column has no numeric value in it. Supplied correct column name?")}
    if(is.null(testChr)) { stop("Chromosome column has no numeric value in it. Supplied correct column name?")}
    if(! is.numeric(testStart)) { stop("Start column has no numeric value in it. Supplied correct column name?")}
    if(! is.numeric(testEnd)) { stop("End column has no numeric value in it. Supplied correct column name?")}
    
    # Add diff column to names we want to keep when merging (comes from function "idSmoothingTargets").
    colNameMerge = c(colNameMerge, "diffs")
    
    registerDoMC(CORES)
    lSmooth = foreach(thisSample = lRaw, .final = function(x) setNames(x, names(lRaw)) ) %dopar% {
        
        thisOut = thisSample
        stillSmoothing = sum(thisOut$smooth)
        while( stillSmoothing > 0 ) {
            # For the while loop:
            # Read lines from thisSample and change in thisOut. Hence for a new iteration I need to sync the two.
            thisSample = thisOut
            
            rleRaw = rle(thisSample$smooth)
            # This takes the indeces of the FALSE chains and adds 1. This should give you the next segment which is TRUE.
            # Two challenges:
            # 1) Last segment always FALSE (see above), hence removal of the last number as this would indicate to a segment outside the df.
            # 2) If it starts with a TRUE segment, this would not be found when looking at the FALSE chains. Hence, adding index 1 manually if chain starts with TRUE.
            indRaw = cumsum(rleRaw$lengths)[ ! rleRaw$values ] + 1
            indRaw = indRaw[ -length(indRaw) ]
            if( rleRaw$values[1] ) { indRaw = c(1, indRaw) }
            
            # loop over start indices of TRUE chains.
            for(i in indRaw) {
                # detect length of segments to smooth. add 1 as the last segment has a FALSE value in it but still belongs to this chain.
                endOfStreak = i + rle(thisSample$smooth[i:nrow(thisSample)])$lengths[1]
                # extract reads
                dfMerge = thisSample[i:endOfStreak,]
                
                # too stupid to make this work with data.table
                newElement = as.data.frame( dfMerge[1,] )
                # Get new end and check first wether valid number.
                newEnd = dfMerge[nrow(dfMerge),][[colNameEnd]]
                if(! is.null(newEnd)) {
                    newElement[[colNameEnd]] = newEnd
                } else {
                    stop("New end coordinate is null. Supplied correct column name?")
                }
                ## Column "segVal" will be dealt with in a minute. Column "diffs" later when running again idSmoothingTargets.
                
                # Merge cn specifically by taking the length of the elements into consideration
                widthWeights = dfMerge[[colNameEnd]] - dfMerge[[colNameStart]]
                newElement[[colNameMerge[1]]] = weighted.mean(dfMerge[[colNameMerge[1]]], widthWeights)
                # Replace all to merge segments with the new merged segment. Later delete duplicated.
                thisOut[i:endOfStreak,] = newElement
            }
            
            # as we have replaced all segments with the new mean segment, we need to remove the duplicates
            thisOut = thisOut[ ! duplicated(thisOut), ]
            # again detect segments which needs smoothing
            thisOut = idSmoothingTargets(thisOut, SMOOTHINGFACTOR, colNameSegVal = colNameMerge[[1]], colNameChr = colNameChr,
                                         IGNOREDELS = IGNOREDELS)
            stillSmoothing = sum(thisOut$smooth)
        }
        
        # after smoothing is finished, change name of cohort
        thisOut$smooth = NULL
        thisOut$diffs = NULL
        return( thisOut )
    }
    
    if( isTRUE(asDf) ) {
        dfSmooth = setDT( rbindlist( lSmooth ) )
        return( dfSmooth )
    } else {
        return( lSmooth )
    }
    
}


## Set everything very close to 2 to 2
dtDat$segVal[ dtDat$segVal > (2-WIGGLE) & dtDat$segVal < (2+WIGGLE) ] = 2
# Merge segments only when two normal follow each other -> SMOOTHINGFACTOR = 0
dtDat = idSmoothingTargets(dtDat, WIGGLE = 0, colNameSegVal = "segVal", colNameChr = "chromosome", IGNOREDELS = IGNOREDELS)
# Split by sample name
lRaw = split(dtDat, dtDat$sample)

## Smooth segments by taking the weighted average of the segVal and their lengths
lSmooth = smoothSegments(lRaw, CORES, SMOOTHINGFACTOR = 0, colNameMerge = "segVal", colNameChr = "chromosome", 
                         colNameStart = "start", colNameEnd = "end", IGNOREDELS = IGNOREDELS, asDf = FALSE)
dtSmooth = rbindlist(lSmooth)

## Write txt and rds
write.table(dtSmooth, paste0(OUTSMOOTHFILE, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
saveRDS(dtSmooth, paste0(OUTSMOOTHFILE, ".rds") )



### PART III: Gather stats about the samples and segments
# Add number of ASCAT_segments and number of CNAs to Kerstin's metadata file
meta=fread(SUMMARYFILE)

# Segments ASCAT detected
dtSegments = data.table(table(dtDat$sample))
meta$segments = dtSegments$N[ match(meta$name, dtSegments$V1) ]

# Segments after smoothing normals
dtSmoothedSegments = data.table(table(dtSmooth$sample))
meta$smoothed_segments = dtSmoothedSegments$N[ match(meta$name, dtSmoothedSegments$V1) ]

# Actual CNAs
dtCNA = dtSmooth[ dtSmooth$segVal != 2, ]
dtCNASegments = data.table(table(dtCNA$sample))
meta$CNAs = dtCNASegments$N[ match(meta$name, dtCNASegments$V1) ]

write.table(meta, SUMMARYFILE, append = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")
