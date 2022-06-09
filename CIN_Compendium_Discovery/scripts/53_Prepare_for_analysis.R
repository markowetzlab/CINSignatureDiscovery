## Intro
args=commandArgs(trailingOnly = TRUE)

MIXMODELSTCGA=args[1]
SIGDEFSTCGA=args[2]
EXPTCGATCGA=args[3]

MIXMODELSOV=args[4]
SIGDEFSOV=args[5]
EXPOVTCGA=args[6]

THRESHOLD=as.numeric(args[7])
NOCN=as.logical(args[8])
RESULTSFILE=args[9]
FUNCTIONS=args[10]

print(paste("MIXMODELSTCGA:", MIXMODELSTCGA))
print(paste("SIGDEFSTCGA:", SIGDEFSTCGA))
print(paste("EXPTCGATCGA:", EXPTCGATCGA))
print(paste("MIXMODELSOV:", MIXMODELSOV))
print(paste("SIGDEFSOV:", SIGDEFSOV))
print(paste("EXPOVTCGA:", EXPOVTCGA))
print(paste("THRESHOLD:", THRESHOLD))
print(paste("NOCN:", NOCN))
print(paste("RESULTSFILE:", RESULTSFILE))
print(paste("FUNCTIONS:", FUNCTIONS))

# ## For testing purposes only
# DATA="/Users/drews01/phd/prjcts/cnsigs2/data"
# CHAPTER1FOLDER="2_OV_signatures_on_TCGA"
# CHAPTER2FOLDER="3_Pancancer_Signatures"
# MIXMODELSTCGA=file.path(DATA, CHAPTER2FOLDER, "2_combined_mixmodels_merged_components.rds")
# SIGDEFSTCGA=file.path(DATA, "5_Signature_Compendium/2_Signature_Compendium_Cosine-0.74.rds")
# # SIGDEFSTCGA=file.path(DATA, CHAPTER2FOLDER, "4_Signatures_KRJ5F9_normalised.rds")
# EXPTCGATCGA=file.path(DATA, "5_Signature_Compendium/3_Exposures_Signature_Compendium_Cosine-0.74.rds")
# # EXPTCGATCGA=file.path(DATA, CHAPTER2FOLDER, "4_Exposures_KRJ5F9_normalised.rds")
# # 
# MIXMODELSOV=file.path(DATA, "Geoffs_data/component_parameters.rds")
# SIGDEFSOV=file.path(DATA, "Geoffs_data/feat_sig_mat.rds")
# EXPOVTCGA=file.path(DATA, CHAPTER1FOLDER, "/Export-matrix_OV_Sigs_on_TCGA-OV_12112019.rds")
# 
# THRESHOLD=0.85
# NOCN=TRUE
# RESULTSFILE="/Users/drews01/phd/prjcts/cnsigs2/results/Signature_Compendium_v5_Cosine-0.74"
# FUNCTIONS="/Users/drews01/phd/prjcts/cnsigs2/finalScripts/53_Prepare_for_analysis_functions.R"


## Packages
# For accessing the OV mixture models
library(flexmix)
# For data handling
library(data.table)
# For cosine function
library(lsa)
# For plotting
library(ggplot2)
library(RColorBrewer)
library(ggthemes)


## Functions
source(FUNCTIONS)

## Rename signatures TBD
expTcga = readRDS(EXPTCGATCGA)
saveRDS(expTcga, paste0(RESULTSFILE, "_Exposures_originalNames.rds"))

newOrder = names(sort(colSums(expTcga > 0.05), decreasing = TRUE))
newNames = paste0(paste0("CS", 1:length(newOrder)), "-", newOrder)

newExpTcga = expTcga[,newOrder]
colnames(newExpTcga) = newNames
saveRDS(newExpTcga, paste0(RESULTSFILE, "_Exposures_bothNames.rds"))

colnames(newExpTcga) = paste0("CS", 1:ncol(newExpTcga))
saveRDS(newExpTcga, paste0(RESULTSFILE, "_Exposures_newNames.rds"))

# Signature definitions
sigTCGA = readRDS(SIGDEFSTCGA)
saveRDS(sigTCGA, paste0(RESULTSFILE, "_Signatures_originalNames.rds"))

if( ! identical(colnames(expTcga), rownames(sigTCGA)) ) { stop("Colnames of exposure file and row names of signature file are not identical!") }

newSigTcga = sigTCGA[newOrder,]
rownames(newSigTcga) = newNames
saveRDS(newSigTcga, paste0(RESULTSFILE, "_Signatures_bothNames.rds"))

rownames(newSigTcga) = paste0("CS", 1:nrow(newSigTcga))
saveRDS(newSigTcga, paste0(RESULTSFILE, "_Signatures_newNames.rds"))


## Get conversion matrix
convMatrix = getConversionMatrix(MODELORIGIN = MIXMODELSTCGA, MODELTARGET = MIXMODELSOV, plotHeatmap = FALSE, uninfPrior = FALSE, noCN = NOCN)

## Liftover and compare signatures
SIGDEFSTCGA=paste0(RESULTSFILE, "_Signatures_newNames.rds")
liftedTCGASigs = liftOverSigs(SIGDEFSTCGA, convMatrix)

# Rows are the 7 OV signatures, columns the TCGA PanCancer signatures
sigCos = compareSigs(liftedTCGASigs, QUERYSIGS = SIGDEFSOV, convMatrix, cosineThresh = THRESHOLD, showPlot = FALSE)

## Compare exposures
EXPTCGATCGA=paste0(RESULTSFILE, "_Exposures_newNames.rds")
expCos = compareExposures(EXPTCGATCGA, EXPOVTCGA, cosineThresh = THRESHOLD, showPlot = FALSE)

## Combine validations
plotBoth = combineInformation(sigCos, expCos, PREFIXCOLS = "OV", PREFIXROWS = "CS", TRESHOLD = 0.85, 
                              CATSIG = "Definition", CATEXP = "Exposure", CATBOTH = "Def+Exp")
    
pdf(paste0(RESULTSFILE, "_compare_to_OV_sigs.pdf"), width = 12, height = 4)
print(plotBoth)
dev.off()

## Convert plot to data.frame for later figure in paper
dtTable = plotBoth$data
dtTable$value = NULL
dtTable$Name = NULL
dtTable$Test = "OV Signatures"

write.table(dtTable, file = paste0(RESULTSFILE, "_table_OV_sigs_comparison.txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)

