# Plot correlations about pancancer and cancerspecific signatures
# Decide on cancerspecific signatures which seem not to be in pancancer set

library(ComplexHeatmap)
library(lsa)
library(data.table)
library(ggthemes)
library(ggplot2)

args=commandArgs(trailingOnly = TRUE)
PANCANCERSIGS=args[1]
print(paste("PANCANCERSIGS:", PANCANCERSIGS))
PANCANCEREXP=args[2]
print(paste("PANCANCEREXP:", PANCANCEREXP))
PATHCANCERSPECIFIC=args[3]
print(paste("PATHCANCERSPECIFIC:", PATHCANCERSPECIFIC))

PATTERNSIGS=args[4]
print(paste("PATTERNSIGS:", PATTERNSIGS))
PATTERNEXP=args[5]
print(paste("PATTERNEXP:", PATTERNEXP))

OUTDIR=args[6]
print(paste("OUTDIR:", OUTDIR))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

COSINETHRESH=as.numeric(args[7])
print(paste("COSINETHRESH:", COSINETHRESH))

CANCERCOLS=args[8]
print(paste("CANCERCOLS:", CANCERCOLS))
SOURCEFUNCTIONS=args[9]
print(paste("SOURCEFUNCTIONS:", SOURCEFUNCTIONS))

## Testing
# PANCANCERSIGS="/Users/drews01/phd/prjcts/cnsigs2/data/3_Pancancer_Signatures/4_Signatures_KRJ5F9_normalised.rds"
# PANCANCEREXP="/Users/drews01/phd/prjcts/cnsigs2/data/3_Pancancer_Signatures/4_Exposures_KRJ5F9_normalised.rds"
# PATHCANCERSPECIFIC="/Users/drews01/phd/prjcts/cnsigs2/data/4_Cancer_specific_Signatures"
# PATTERNSIGS="4_Signatures"
# PATTERNEXP="4_Exposures"
# OUTDIR="/Users/drews01/phd/prjcts/cnsigs2/data/5_Signature_Compendium"
# dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
# COSINETHRESH=0.74
# CANCERCOLS="/Users/drews01/phd/prjcts/cnsigs2/data/metadata/TCGA_colour_scheme_Lydia.txt"
# SOURCEFUNCTIONS="/Users/drews01/phd/prjcts/cnsigs2/finalScripts/51_Combine_PC_and_CS_Signatures_functions.R"


## Input from hard drive
source(SOURCEFUNCTIONS)


# Load pan-cancer data
pcSigs = readSigs(PANCANCERSIGS, name = "TCGA")
pcExp = t(readExp(PANCANCEREXP, name = "TCGA"))


# Cancer-specific signatures
allCSigFiles = list.files( path = PATHCANCERSPECIFIC, pattern = PATTERNSIGS, full.names = TRUE,recursive = TRUE)
allCSigFiles = allCSigFiles[ grep( "normalised", allCSigFiles, invert = TRUE ) ]
allCSigFiles = allCSigFiles[ grep( "not_used", allCSigFiles, invert = TRUE ) ]
allCSigFiles = allCSigFiles[ grep( "summary", allCSigFiles, invert = TRUE ) ]
lCSigs = lapply(allCSigFiles, function(thisFile) {
    
    thisName = basename(dirname(thisFile))
    thisSig = readSigs(thisFile, name = thisName)
    return(thisSig)
    
} )
dfCS = do.call(rbind, lCSigs)
matAll = as.matrix( rbind(dfCS, pcSigs) )


# Cancer-specific exposures
allCExpFiles = list.files( path = PATHCANCERSPECIFIC, pattern = PATTERNEXP, full.names = TRUE,recursive = TRUE)
allCExpFiles = allCExpFiles[ grep( "normalised", allCExpFiles ) ]
allCExpFiles = allCExpFiles[ grep( "not_used", allCExpFiles, invert = TRUE ) ]
lCExp = lapply(allCExpFiles, function(thisFile) {
    
    thisName = basename(dirname(thisFile))
    thisExp = t(readSigs(thisFile, name = thisName))
    return(thisExp)
    
} )
names(lCExp) = basename(dirname(allCExpFiles))


# Prepare colours
cancerCols = read.table(CANCERCOLS, header = FALSE, comment.char = "|")
vCols = as.character(cancerCols$V2)
names(vCols) = cancerCols$V1
vCols = c(vCols, "TCGA" = "#636363")



## Plot hairball and hierarchical clustering
outfileHairball = file.path(OUTDIR, paste0("0_Hairball_PC_and_CS_sigs_cosineSim"))
plotHairball(matAll, outfileHairball, numPCSigs = nrow(pcSigs), corThreshold = COSINETHRESH, 
             cancerCols = cancerCols, corMethod = "cosine", plot2File = TRUE)

corMat = cosine(t(matAll))
outfileHierclust = file.path(OUTDIR, paste0("0_Hierarchical_clustering_PC_and_CS_Sigs_cosineSim.pdf"))
pdf(outfileHierclust, width = 12, height = 12)
Heatmap(corMat, show_row_names = FALSE, show_column_names = FALSE, name = "Cosine_Sim")
dev.off()

corMatFilt = corMat
corMatFilt[ corMatFilt < COSINETHRESH ] = 0
outfileHierclust2 = file.path(OUTDIR, paste0("0_Hierarchical_clustering_PC_and_CS_Sigs_cosineSim_filt-", COSINETHRESH,".pdf"))
pdf(outfileHierclust2, width = 12, height = 12)
Heatmap(corMatFilt, show_row_names = FALSE, show_column_names = FALSE, name = "Cosine_Sim")
dev.off()



## Step 1: Number of exposed samples per signature
dtExposedSamples = samplesPerSignature(lCExp, EXPTRESH = 0.05)

outfileHistSigExp = file.path(OUTDIR, "0_Histogram_CS_and_PC_Signature_Exposures.pdf")
pdf(outfileHistSigExp, width = 12, height = 8)
print( ggplot(dtExposedSamples, aes(x = Samples, fill = Cancer)) + geom_histogram() + scale_x_log10() + 
    scale_fill_manual(values = vCols) + theme_tufte(base_family = "", base_size = 24) + ylab("Count") + 
    ggtitle("Frequency of signature exposure") + guides(fill = guide_legend(ncol= 1)) )
dev.off()


## Step 2: CS sig vs PC sig
filtMatAll = removePCSimilarSigs(matAll, COSINETHRESH = COSINETHRESH)

## Step 3: CS sigs vs CS sig (retain signature with most exposed samples)
csSigs = removeCSSimilarSigs(filtMatAll, dtExposedSamples, COSINETHRESH = COSINETHRESH)

## Step 4: Merge PC and surviving CS sigs to create Pan-cancer compendium signatures (PCCS)
matPCSC = rbind(pcSigs, csSigs)

## Step 4: NNLS (not needed since no two PC signatures are combining to a CS sig)
# library(nnls)
# ALLOWEDERROR=0.1
# for(index in 1:nrow(csSigs)) {
#     
#     thisSig = csSigs[ index, ]
#     thisSigName = rownames(csSigs)[index]
#     print(thisSigName)
#     
#     LIMIT = sum(thisSig)*ALLOWEDERROR
#     
#     # Try to fit pcSigs together to form the CS sig
#     # Assumption: Pancancer signatures have more data in them and thus are theoretically able to detect more detailed processes
#     # Assumption2: CS signatures may not be as fine-grained as PC sigs and thus may be puzzled together by PC sigs.
#     model=nnls(t(pcSigs), thisSig)
#     sumRes = sum(abs(model$residuals))
#     print(sumRes)
#     print(sum(model$x != 0))
# }


## Plot the new signatures
outfileHeatmaps = file.path(OUTDIR, paste0("1_Heatmaps_Cosine-", COSINETHRESH, ".pdf"))
pdf(outfileHeatmaps, width = 14, height = 12)
Heatmap(matPCSC, cluster_columns = FALSE, cluster_rows = FALSE, column_title = "Pan-cancer signature compendium")

mOther = matPCSC
mOther = apply(mOther, 2, function(y) y/sum(y))
Heatmap(mOther, cluster_columns = FALSE, column_title = "PCSC (normalised per component)")

Heatmap(pcSigs, cluster_columns = FALSE, cluster_rows = FALSE, column_title = "Pan-cancer signatures")
Heatmap(csSigs, cluster_columns = FALSE, cluster_rows = FALSE, column_title = "Cancer-specific signatures")
dev.off()

# outfileExpl = file.path(OUTDIR, paste0("2_Signature_Compendium_Cosine-", COSINETHRESH, "_explanations"))
# write.table(dtExplain, paste0(outfileExpl, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
# saveRDS(dtExplain, paste0(outfileExpl, ".rds"))

outfileSigs = file.path(OUTDIR, paste0("2_Signature_Compendium_Cosine-", COSINETHRESH, ".rds"))
saveRDS(matPCSC, outfileSigs)

outCSS = file.path(OUTDIR, paste0("2_Signatures_CSS_Cosine-", COSINETHRESH, ".rds"))
saveRDS(csSigs, outCSS)

outPCS = file.path(OUTDIR, paste0("2_Signatures_PCS_Cosine-", COSINETHRESH, ".rds"))
saveRDS(pcSigs, outPCS)


#### TODO
# ### Final plot: Alluvial plot of how signatures are classified
# ## Combine information for alluvial
# dtExposedSamples$PCCS = "Rejected"
# dtExposedSamples$PCCS[ rownames(dtExposedSamples) %in% rownames(matPCCS) ] = "Admitted"
# dtAlluv = data.table(table(dtExposedSamples$Cancer, dtExposedSamples$PCCS))
# dtAlluv = dtAlluv[order(dtAlluv$N),]
# ggplot(dtAlluv,
#        aes(y = N, axis1 = V1, axis2 = V2)) +
#     geom_alluvium(aes(fill = V2), width = 1/12) +
#     geom_stratum(width = 1/12, fill = "black", color = "grey")
# 
# dtAlluv = data.table(table(dtExposedSamples$Cancer))
# dtAlluv$Strata = "Before"
# 
# 
# library(ggalluvial)
# 
# ggplot(as.data.frame(UCBAdmissions),
#        aes(y = Freq, axis1 = Gender, axis2 = Dept)) +
#     geom_alluvium(aes(fill = Admit), width = 1/12) +
#     geom_stratum(width = 1/12, fill = "black", color = "grey") +
#     geom_label(stat = "stratum", infer.label = TRUE) +
#     scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
#     scale_fill_brewer(type = "qual", palette = "Set1") +
#     ggtitle("UC Berkeley admissions and rejections, by sex and department")



