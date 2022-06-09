#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 10) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

## General configs
# General
PATHMAIN=args[1]
print(paste("PATHMAIN:", PATHMAIN))
INPUTECNF=args[2]
print(paste("INPUTECNF:", INPUTECNF))
OUTPATH=args[3]
dir.create(OUTPATH, showWarnings = FALSE, recursive = TRUE)
print(paste("OUTPATH:", OUTPATH))
FEATURE=args[4]
print(paste("FEATURE:", FEATURE))

# Mix model specific
FITDIST=args[5]
print(paste("FITDIST:", FITDIST))
MINCOMP=as.numeric( args[6] )
print(paste("MINCOMP:", MINCOMP))
MAXCOMP=as.numeric( args[7] )
print(paste("MAXCOMP:", MAXCOMP))
PRIOR=as.numeric( args[8] )
print(paste("PRIOR:", PRIOR))
NITER=as.numeric( args[9] )
print(paste("NITER:", NITER))
DECISION=args[10]
print(paste("DECISION:", DECISION))
SEED=as.numeric(args[11])
print(paste("SEED:", SEED))

# Set seed
if( is.na(SEED) ) {
    SEED = 77777
}
print(paste("SEED:", SEED))


# ## For testing
# PATHMAIN="/home/drews01/data/phd/cnsigs2/data/main_functions.R"
# INPUTECNF="/home/drews01/data/phd/cnsigs2/data/3_Pancancer_Signatures/1_tcga_filtered_ecnf.rds"
# OUTPATH="/home/drews01/data/phd/cnsigs2/data/3_Pancancer_Signatures/"
# # OUTPATH="/home/drews01/phd/prjcts/cnsigs2/data/3_TCGA_specific_signatures"
# FEATURE="bp10MB"
# FITDIST="pois"
# MINCOMP=2
# MAXCOMP=5
# PRIOR=0.001
# NITER=1000
# DECISION="BIC"
# SEED=77777

source(PATHMAIN)
# To be sure this is properly loaded
# Poisson
library(flexmix)
# Gaussian
library(mclust)
library(data.table)
library(ggplot2)
library(cowplot)

# Load data
tcga.ecnf = readRDS(INPUTECNF)

# Remove old files
FILENAME = paste0( "2_fit_", FEATURE, "_prior-", PRIOR, "_min-", MINCOMP, "_max-", MAXCOMP, "_modeldec-", DECISION )
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
tryCatch(
    {
    # most of the options only needed for Poisson mixture models (as run with flexmix).
    # Gaussian mixture models done with Mclust.
    fittedDist = fitComponent2( dat, dist = FITDIST, min_prior = PRIOR, min_comp = MINCOMP, 
                                   max_comp = MAXCOMP, niter = NITER, model_selection = DECISION, seed = SEED )
    saveRDS( fittedDist,  paste0( OUTNAME, ".RDS") )

    },
    error = function(cond) {
        file.remove( startFlag )
        write.table( dummy, file = paste0( OUTNAME, ".FAILED" ), col.names = FALSE )

        message( "This run failed with the following error message:" )
        message( cond )
    },
    warning = function(cond) {
        file.remove( startFlag )
        write.table( dummy, file = paste0( OUTNAME, ".WARNING-seelogs" ), col.names = FALSE )

        message( "This run had a warning message:" )
        message( cond )
    }
)

# Set flags
# Save output
if( ! exists("fittedDist") ) {
    file.remove( startFlag )
    write.table( dummy, file = paste0( OUTNAME, ".FAILED" ), col.names = FALSE )
    
    stop("Fitting the mixture models failed. The script will stop here.")
} else {
    # Set the right flags and save output
    file.remove( startFlag )
    write.table( dummy, file = paste0( OUTNAME, ".SUCCESS" ), col.names = FALSE )
}

# Convert cluster components to a table for easier validation
if( FITDIST == "pois" ) {
    # Poisson distribution -> flexmix object
    listComp = unlist(fittedDist@components)
    allComps = sapply( listComp, function(thisComp) { unlist( thisComp@parameters ) })
    dtComp = data.table( cluster = seq(1:length(allComps)), lambda = sort(allComps))
} else {
    # Normal distribution -> Mclust object
    dtComp = data.table( cluster = names(fittedDist$parameters$mean), mean = fittedDist$parameters$mean,
                         sd = sqrt(fittedDist$parameters$variance$sigmasq) )
}

write.table( dtComp, file = paste0( OUTNAME, ".txt"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )
