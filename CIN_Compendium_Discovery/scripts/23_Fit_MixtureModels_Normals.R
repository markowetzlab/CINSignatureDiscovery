#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

## General configs
# General
INPUTECNF=args[1]
print(paste("INPUTECNF:", INPUTECNF))
OUTPATH=args[2]
dir.create(OUTPATH, showWarnings = FALSE, recursive = TRUE)
print(paste("OUTPATH:", OUTPATH))
FEATURE=args[3]
print(paste("FEATURE:", FEATURE))

# Mix model specific
MINCOMP=as.numeric( args[4] )
print(paste("MINCOMP:", MINCOMP))
MAXCOMP=as.numeric( args[5] )
print(paste("MAXCOMP:", MAXCOMP))
RANGE=MINCOMP:MAXCOMP
ITER=as.numeric( args[6] )
print(paste("ITER:", ITER))
CORES=as.numeric( args[7] )
print(paste("CORES:", CORES))


# ## For testing
# PATHMAIN="/home/drews01/data/phd/cnsigs2/data/main_functions.R"
# INPUTECNF="/home/drews01/data/phd/cnsigs2/data/3_TCGA_specific_signatures/1_tcga_filtered_ecnf.rds"
# OUTPATH="/home/drews01/data/phd/cnsigs2/data/3_TCGA_specific_signatures"
# # OUTPATH="/home/drews01/phd/prjcts/cnsigs2/data/3_TCGA_specific_signatures"
# FEATURE="segsize"
# MINCOMP=1
# MAXCOMP=50
# RANGE=MINCOMP:MAXCOMP
# ITER=5
# CORES=5

library(mclust)
library(data.table)
library(foreach)
library(doMC)


# Functions
calcEntau = function(thisModel) {
    matZ = thisModel$z
    logMatZ = log(matZ)
    logMatZ[ is.infinite(logMatZ) ] = -1000
    entau = -sum( matZ*logMatZ )
    return(entau)
}

iclbic = function(thisModel, entau) {
    iclbic = -2*thisModel$loglik + 2*entau + thisModel$df*log(thisModel$n)
    return(iclbic)
}


# Load data
tcga.ecnf = readRDS(INPUTECNF)

# Remove old files
FILENAME = paste0( "2_fit_", FEATURE, "_iter-", ITER, "_min-", MINCOMP, "_max-", MAXCOMP )
OLDFILES = list.files(path = OUTPATH, pattern = paste0(FILENAME, "*") )
if( length(OLDFILES) > 0 ) {
    for(THISFILE in OLDFILES) {
        file.remove( file.path(OUTPATH, THISFILE) )        
    }
}

OUTNAME = file.path( OUTPATH, FILENAME )

# Fit model. Take time and set flag to see whether job dies.
dummy = data.frame()
startFlag = paste0( OUTNAME, ".STARTED" )
write.table( dummy, file = startFlag, col.names = FALSE )

dat = as.numeric( tcga.ecnf[[ FEATURE ]][[2]] )

registerDoMC(cores = CORES)
lResultsAllModels = foreach(thisIter=1:ITER) %dopar% {
    
    # debugging
    print(paste("This iteration:", thisIter))
    
    lFullOut = lapply(RANGE, function(thisK) {
        
        print(paste("  This K:", thisK))
        
        # Hashtag reproducibility
        thisSeed = ITER*thisIter + thisK
        print(paste("    This seed:", thisSeed))
        set.seed(thisSeed)
        
        # Fit and extract most needed information. I cannot store the model itself as it would consume too much disk space.
        thisFit = Mclust(dat, modelNames = "V", G = thisK)
        # In rare instances it doesn't find a solution and thus I run it again with the inferior model where the variances are the same.
        if( class(thisFit) != "Mclust" ) {
		    set.seed(thisSeed)
		    thisFit = Mclust(dat, G = thisK)
        }
        
	    thisEN = calcEntau(thisFit)
        thisOut = c(thisK, thisFit$loglik, thisFit$df, thisFit$n, thisEN, abs(as.numeric(thisFit$bic)), 
                    abs(as.numeric(icl(thisFit))), iclbic(thisFit, thisEN) )
        
        # For controlling the max memory on the cluster.
        rm(thisFit)
        invisible(gc())
        
        return( thisOut )
        
    } )
    
    # Prepare output
    dtFullOut = data.table( do.call("rbind", lFullOut) )
    sampleID = paste0(FEATURE, "-", thisIter)
    dtFullOut$sampleID = sampleID
    names(dtFullOut) = c("fitK", "logLH", "df", "n", "EN", "BIC", "ICL", "BICICL", "sampleID" )
    
    # Get best K for BIC, ICL and BIC-ICL for storing in simple output
    dtSimpleOut = c(dtFullOut$fitK[ which.min(dtFullOut$BIC) ], dtFullOut$fitK[ which.min(dtFullOut$ICL) ], 
                    dtFullOut$fitK[ which.min(dtFullOut$BICICL) ], sampleID)
    names(dtSimpleOut) = c("BIC", "ICL", "BICICL", "sampleID")
    
    # Return output
    lOut = list()
    lOut$full = dtFullOut
    lOut$simple = dtSimpleOut
    
    return( lOut )
    
}

print("Modelling finished. Now converting to output and saving.")

# Combine simples and fulls together to store independently
lResultsFull = lapply(lResultsAllModels, function(thisResult) {
    return( thisResult$full )
})

dtResultsFull = rbindlist(lResultsFull)

lResultsSimple = lapply(lResultsAllModels, function(thisResult) {
    return( thisResult$simple )
})

dtResultsSimple = as.data.table( do.call("rbind", lResultsSimple) )

OUTNAMEFULL = paste0(OUTNAME, "_full.RDS")
OUTNAMESIMPLE = paste0(OUTNAME, "_simple.RDS")

saveRDS(dtResultsFull, OUTNAMEFULL)
saveRDS(dtResultsSimple, OUTNAMESIMPLE)


OUTNAMESIMPLE = paste0(OUTNAME, "_simple.txt")

write.table(dtResultsSimple, OUTNAMESIMPLE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

file.remove( startFlag )

print("Finished.")
