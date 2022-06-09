## Contrary to its name this script also looks at the differences the new feature extraction methods has on the OV exposures.

## There are large differences of OV signature exposures on TCGA-OV when I use Geoff's CN feature extraction method or mine. The differences
## are exacerbated by using older (pre-optimised ASCAT) TCGA segmentation.
## This script aims to find out where this discrepancy is coming from. So it goes step by step through the different combinations of data and
## CN extraction methods.
##
## There are three data sets:
### Old segmentation (pre-optimised ASCAT) 
### New optimised ASCAT
### New optimised ASCAT + smoothing of normal segments (collapsing and merging)
##
## There are three types of CN feature extraction:
### Old method (Geoff's Nat Gen)
### New method (Re-writing by me for large-cohort optimisation and removal of changepoint bug [ignoring first segment])
### New method + remove normals
##
## I don't have the raw data for the old method just Geoff's extracted CN features, therefor I cannot do all permutations of comparing
## signature exposures. But I will start with old/old and new/old. Next one is new/old and new/new.

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 12) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

## Preface
library(data.table)
# Melt
library(reshape)
# Plotting
library(ggplot2)
library(ggthemes)
library(ComplexHeatmap)

PREPATH=args[1]
print(paste("PREPATH:", PREPATH))
OUT=args[2]
print(paste("OUT:", OUT))
OVSIGS=args[3]
print(paste("OVSIGS:", OVSIGS))
GEOFFPARAMS=args[4]
print(paste("GEOFFPARAMS:", GEOFFPARAMS))
CORES=as.numeric(args[5])
print(paste("CORES:", CORES))
OLDMETHODS=args[6]
print(paste("OLDMETHODS:", OLDMETHODS))
NEWMETHODS=args[7]
print(paste("NEWMETHODS:", NEWMETHODS))
META=args[8]
print(paste("META:", META))
OLDCN=args[9]
print(paste("OLDCN:", OLDCN))
NEWRAW=args[10]
print(paste("NEWRAW:", NEWRAW))
NEWPLUS=args[11]
print(paste("NEWPLUS:", NEWPLUS))
DATE=args[12]
print(paste("DATE:", DATE))

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ## Testing purposes
BASE="/Users/drews01/phd/prjcts/cnsigs2/data"
# BASE="/Users/drews01/data/phd/cnsigs2/data"
PREPATH="/Users/drews01/phd/prjcts/cnsigs2/data/Geoffs_data/"
OUT=file.path(BASE, "2_OV_signatures_on_TCGA_temp")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
OVSIGS=file.path(BASE, "Geoffs_data/feat_sig_mat.rds")
GEOFFPARAMS=file.path(BASE, "Geoffs_data/component_parameters.rds")
CORES=2
OLDMETHODS=file.path(BASE, "../src/cnsignatures_vanilla/main_functions.R")
NEWMETHODS=file.path(BASE, "main_functions.R")
META=file.path(BASE, "metadata/summary.ascatTCGA.penalty70.txt")
OLDCN=file.path(BASE, "validationdata/data_Geoff/data/tcga_CN_features.rds")
NEWRAW=file.path(BASE, "rawdata/combined.ascat.segments.filt.rds")
NEWPLUS=file.path(BASE, "rawdata/combined.ascat.segments.smoothednormals.rds")
DATE="12112019"


ovSigs=readRDS(OVSIGS)
gm_params=readRDS(GEOFFPARAMS)
meta = fread(META)

ovSamps = substr(meta$name[ meta$cancer_type == "OV" & ! is.na(meta$CNAs) ],1,12)

# Function
filterAndOrder = function(set1, set2) {
    set1 = set1[ rownames(set1) %in% rownames(set2), ]
    set2 = set2[ rownames(set2) %in% rownames(set1), ]
    set2 = set2[ match(rownames(set1), rownames(set2)), ]
    # Check if they are now identical
    if(! identical(rownames(set1), rownames(set1))) stop("Sets cannot be matched.")
    return(list(set1 = set1, set2 = set2))    
}

compareExposures = function(set1, set2, nameSet1, nameSet2) {
    
    # Filter and match order
    lSets = filterAndOrder(set1, set2)
    
    # Combine info
    m1 = melt(lSets$set1)
    m2 = melt(lSets$set2)
    if(! identical(m1$X1, m2$X1) & identical(m1$X2, m2$X2)) stop("Melting resulted in error.")
    m1$NewOld = m2$value
    colnames(m1) = c("Samples", "Signature", nameSet1, nameSet2)
    
    returnPlot = ggplot(m1, aes_string(x = nameSet1, y = nameSet2, colour = "Signature")) + geom_point(size = 0.75, alpha = 0.6) + 
        facet_wrap(.~Signature, scales = "free_x") + geom_abline(intercept = 0, slope = 1) + theme_tufte(base_family = "") + 
        ggtitle("OV sigs on TCGA-OV from two CN feature extraction methods")
    return(returnPlot)
}

compareExposuresHM = function(set1, set2, nameSet1, nameSet2) {
    
    # Filter and match order
    lSets = filterAndOrder(set1, set2)
    
    lPlots = list()
    lPlots$set1 = Heatmap(lSets$set1, show_row_names = FALSE, cluster_rows = TRUE, 
                          cluster_columns = FALSE, column_title = paste(nameSet1))
    lPlots$set2 = Heatmap(lSets$set2, show_row_names = FALSE, cluster_rows = TRUE, 
                          cluster_columns = FALSE, column_title = paste(nameSet2))
    lPlots$diff = Heatmap(lSets$set1 - lSets$set2, show_row_names = FALSE, cluster_rows = TRUE, 
                          cluster_columns = FALSE, column_title = paste(nameSet1, "minus", nameSet2))
    
    return(lPlots)
}



## Comparison I: old/old vs new/old
# old/old
source(OLDMETHODS)
oldCN = readRDS(OLDCN)
oldSxC = generateSampleByComponentMatrix(oldCN, all_components = gm_params, cores = CORES)
oldSxC = oldSxC[ ,match(rownames(ovSigs), colnames(oldSxC)) ]
oldOldExp = t( quantifySignatures(oldSxC, component_by_signature = ovSigs) )
rownames(oldOldExp) = substr(rownames(oldOldExp), 1, 12)

# new/old
# Load raw TCGA ASCAT data and get OV samples
newRaw = readRDS(NEWRAW)
newRaw = newRaw[ newRaw$sample %in% ovSamps, ]
newRaw$sample = factor(newRaw$sample)
# ECNF and SxC matrix derivation
newRaw = as.data.frame(newRaw)
newRawList = split(newRaw, newRaw$sample)
newOldCN = extractCopynumberFeatures(newRawList, cores = CORES)
newOldSxC = generateSampleByComponentMatrix(newOldCN, all_components = gm_params, cores = CORES)
newOldSxC = newOldSxC[ ,match(rownames(ovSigs), colnames(newOldSxC)) ]
newOldExp = t( quantifySignatures(newOldSxC, component_by_signature = ovSigs) )
rownames(newOldExp) = substr(rownames(newOldExp), 1, 12)

# Compare both
pOOvNO = compareExposures(set1 = oldOldExp, set2 = newOldExp, nameSet1 = "OldDatOldMethod", nameSet2 = "NewDatOldMethod")
lHeatmaps = compareExposuresHM(set1 = oldOldExp, set2 = newOldExp, nameSet1 = "OldDatOldMethod", nameSet2 = "NewDatOldMethod")
 


## Comparison II: old/old vs new/new
source(NEWMETHODS)
# new/new
# ECNF and SxC matrix derivation
newNewCN = extractCopynumberFeatures(newRawList, cores = CORES, prePath=PREPATH, rmNorm = FALSE)
newNewSxC = generateSampleByComponentMatrix(newNewCN, all_components = gm_params, cores = CORES)
newNewSxC = newNewSxC[ ,match(rownames(ovSigs), colnames(newNewSxC)) ]
newNewExp = t( quantifySignatures(newNewSxC, component_by_signature = ovSigs) )
rownames(newNewExp) = substr(rownames(newNewExp), 1, 12)

# Compare both
# Sanity check
pNNvNO = compareExposures(set1 = newNewExp, set2 = newOldExp, nameSet1 = "NewDatNewMethod", nameSet2 = "NewDatOldMethod")
# Compare to old
pOOvNN = compareExposures(set1 = oldOldExp, set2 = newNewExp, nameSet1 = "OldDatOldMethod", nameSet2 = "NewDatNewMethod")



## Comparison III: old/old vs new+/old
source(OLDMETHODS)
# Load raw TCGA ASCAT data and get OV samples
newPlusRaw = readRDS(NEWPLUS)
newPlusRaw = newPlusRaw[ newPlusRaw$sample %in% ovSamps, ]
newPlusRaw$sample = factor(newPlusRaw$sample)
# ECNF and SxC matrix derivation
newPlusRaw = as.data.frame(newPlusRaw)
newPlusRawList = split(newPlusRaw, newPlusRaw$sample)
newPlusOldCN = extractCopynumberFeatures(newPlusRawList, cores = CORES)
newPlusOldSxC = generateSampleByComponentMatrix(newPlusOldCN, all_components = gm_params, cores = CORES)
newPlusOldSxC = newPlusOldSxC[ ,match(rownames(ovSigs), colnames(newPlusOldSxC)) ]
newPlusOldExp = t( quantifySignatures(newPlusOldSxC, component_by_signature = ovSigs) )
rownames(newPlusOldSxC) = substr(rownames(newPlusOldSxC), 1, 12)

# Compare both
pOOvPO = compareExposures(set1 = oldOldExp, set2 = newPlusOldExp, nameSet1 = "OldDatOldMethod", nameSet2 = "NewPlusDatOldMethod")
pPOvNO = compareExposures(set1 = newPlusOldExp, set2 = newOldExp, nameSet1 = "NewPlusDatOldMethod", nameSet2 = "NewDatOldMethod")
pPOvNN = compareExposures(set1 = newPlusOldExp, set2 = newNewExp, nameSet1 = "NewPlusDatOldMethod", nameSet2 = "NewDatNewMethod")


## Comparison IV: old/old vs new+/new+
source(NEWMETHODS)
# new/new
# ECNF and SxC matrix derivation
newPlusNewPlusCN = extractCopynumberFeatures(newPlusRawList, cores = CORES, prePath=PREPATH, rmNorm = TRUE)
newPlusNewPlusSxC = generateSampleByComponentMatrix(newPlusNewPlusCN, all_components = gm_params, cores = CORES)
newPlusNewPlusSxC = newPlusNewPlusSxC[ ,match(rownames(ovSigs), colnames(newPlusNewPlusSxC)) ]
newPlusNewPlusExp = t( quantifySignatures(newPlusNewPlusSxC, component_by_signature = ovSigs) )
rownames(newPlusNewPlusExp) = substr(rownames(newPlusNewPlusExp), 1, 12)

# Compare both
pOOvPP = compareExposures(set1 = oldOldExp, set2 = newPlusNewPlusExp, nameSet1 = "OldDatOldMethod", nameSet2 = "NewPlusDatNewPlusMethod")

pPPvPO = compareExposures(set1 = newPlusNewPlusExp, set2 = newPlusOldExp, 
                          nameSet1 = "NewPlusDatNewPlusMethod", nameSet2 = "NewPlusDatOldMethod")
pPPvNO = compareExposures(set1 = newPlusNewPlusExp, set2 = newOldExp, nameSet1 = "NewPlusDatNewPlusMethod", nameSet2 = "NewDatOldMethod")
pPPvNN = compareExposures(set1 = newPlusNewPlusExp, set2 = newNewExp, nameSet1 = "NewPlusDatNewPlusMethod", nameSet2 = "NewDatNewMethod")



## Join all data.frames and plot old/old versus all other versions
# Currently have old/old, new/old, new/new, new+/old, new+/new+
dtSummary = data.table(pOOvNO$data)
temp = data.table(pOOvNN$data)
identical(as.character(dtSummary$Samples), as.character(temp$Samples))
dtSummary$NewDatNewMethod = temp$NewDatNewMethod

temp = data.table(pOOvPO$data)
identical(as.character(dtSummary$Samples), as.character(temp$Samples))
dtSummary$NewPlusDatOldMethod = temp$NewPlusDatOldMethod

temp = data.table(pOOvPP$data)
identical(as.character(dtSummary$Samples), as.character(temp$Samples))
dtSummary$NewPlusDatNewPlusMethod = temp$NewPlusDatNewPlusMethod

## Plot scatter plots side by side versus the old/old one.
mSummary = melt(dtSummary, id.vars = c("Samples", "Signature", "OldDatOldMethod"))

pSummary = ggplot(mSummary, aes(x = value, y = OldDatOldMethod, colour = Signature)) + geom_point(size = 0.75, alpha = 0.6) + 
    facet_grid(variable ~ Signature, scales = "free_x") + geom_abline(intercept = 0, slope = 1) + theme_tufte(base_family = "") + 
    ggtitle("OV sigs on TCGA-OV from multiple CN feature extraction methods") + ylab("Exposure of old data and old method") +
    xlab("Exposure of other data and methods")


## Do the same with correlation scores.
lCors = lapply(unique(dtSummary$Signature), function(sig) {
    
    thisSig = dtSummary[ dtSummary$Signature == sig, ]
    mCor = cor(thisSig[,c(-1,-2)], method = "pearson")
    mCor[ lower.tri(mCor, diag = TRUE) ] = NA
    dtCor = melt(mCor)
    dtCor = dtCor[ ! is.na(dtCor$value), ]
    dtCor$Signature = sig
    return(dtCor)
} )
dtCors = do.call(rbind, lCors)

pCors = ggplot(dtCors, aes(x = Signature, y = value, colour = Signature)) + geom_point(size = 4, alpha = 0.6) + 
    theme_tufte(base_family = "") + ylab("Pearson Correlation") + 
    ggtitle("Correlations between multiple CN feature extraction methods and two data sets")

## Let's look at the individual exposure variances.
dtSummary$SDs = apply(dtSummary[,c(-1,-2)], 1, sd)
pSDs = ggplot(dtSummary, aes(x = Signature, y = SDs)) + geom_jitter(size = 2, alpha = 0.6, width = 0.2, height = 0) + 
    geom_boxplot(outlier.colour = NA) + theme_tufte(base_family = "") + ylab("Standard deviations of exposures") + 
    ggtitle("Sample exposures are relatively stable across multiple CN feature extraction methods")

## Save all output
pdf(file.path(OUT, "0_Summary_Scatterplots.pdf"), width = 12, height = 8); print(pSummary); dev.off()
saveRDS(dtSummary, file.path(OUT, "0_data_Summary.rds"))
pdf(file.path(OUT, "1_Summary_Correlations_and_SDs.pdf"), width = 8, height = 8); print(pCors); print(pSDs); dev.off()

pdf(file.path(OUT, "2_Detail_Old_versus_Other.pdf"), width = 8, height = 8)
print(pOOvNO); print(pOOvNN); print(pOOvPO); print(pOOvPP); dev.off()

pdf(file.path(OUT, "3_Detail_Other_versus_Other.pdf"), width = 8, height = 8)
print(pNNvNO); print(pPOvNO); print(pPOvNN); print(pPPvNO); print(pPPvNN); print(pPPvPO); dev.off()



#### This is the actual output => OV sig exposures on TCGA with new methods
## Save signature exposures of new+/new+
# For Sarwah save the input matrix
saveRDS(newPlusNewPlusSxC, file.path(OUT, paste0("Export-input_matrix_of_TCGA-OV_", DATE, ".rds")))
write.table(newPlusNewPlusSxC, file.path(OUT, paste0("Export-input_matrix_of_TCGA-OV_", DATE, ".txt")), quote = FALSE, sep = "\t")

# To get the full names again
newPlusNewPlusExp = t( quantifySignatures(newPlusNewPlusSxC, component_by_signature = ovSigs) )
saveRDS(newPlusNewPlusExp, file.path(OUT, paste0("Export-matrix_OV_Sigs_on_TCGA-OV_", DATE, ".rds")))
write.table(newPlusNewPlusExp, file.path(OUT, paste0("Export-matrix_OV_Sigs_on_TCGA-OV_", DATE, ".txt")), quote = FALSE, sep = "\t")

mNPNPExp = melt(newPlusNewPlusExp)
saveRDS(mNPNPExp, file.path(OUT, paste0("Export-df_OV_Sigs_on_TCGA-OV_", DATE, ".rds")))
write.table(mNPNPExp, file.path(OUT, paste0("Export-df_OV_Sigs_on_TCGA-OV_", DATE, ".txt")), quote = FALSE, sep = "\t", row.names = FALSE)

# Get absolute exposures needed for dCIN detection.
signature_by_sample = t(YAPSA::LCD(t(newPlusNewPlusSxC), 
                                   YAPSA:::normalize_df_per_dim(ovSigs, 2)))
saveRDS(signature_by_sample, file.path(OUT, paste0("Export-matrix_OV_Sigs_on_TCGA-OV_rawExposures_", DATE, ".rds")))


## Derive exposures on all TCGA.
source(NEWMETHODS)
newPlusRaw = readRDS(NEWPLUS)
# ECNF and SxC matrix derivation
newPlusRaw = as.data.frame(newPlusRaw)
newPlusRaw$sample = factor(newPlusRaw$sample)
newPlusRawList = split(newPlusRaw, newPlusRaw$sample)
newPlusNewPlusCN = extractCopynumberFeatures(newPlusRawList, cores = CORES, prePath=PREPATH, rmNorm = TRUE)
newPlusNewPlusSxC = generateSampleByComponentMatrix(newPlusNewPlusCN, all_components = gm_params, cores = CORES)
newPlusNewPlusSxC = newPlusNewPlusSxC[ ,match(rownames(ovSigs), colnames(newPlusNewPlusSxC)) ]
newPlusNewPlusExp = t( quantifySignatures(newPlusNewPlusSxC, component_by_signature = ovSigs) )

saveRDS(newPlusNewPlusExp, file.path(OUT, paste0("Export-matrix_OV_Sigs_on_TCGA_", DATE, ".rds")))
write.table(newPlusNewPlusExp, file.path(OUT, paste0("Export-matrix_OV_Sigs_on_TCGA_", DATE, ".txt")), quote = FALSE, sep = "\t", row.names = FALSE)

mNPNPExp = melt(newPlusNewPlusExp)
saveRDS(mNPNPExp, file.path(OUT, paste0("Export-df_OV_Sigs_on_TCGA_", DATE, ".rds")))
write.table(mNPNPExp, file.path(OUT, paste0("Export-df_OV_Sigs_on_TCGA_", DATE, ".txt")), quote = FALSE, sep = "\t", row.names = FALSE)
