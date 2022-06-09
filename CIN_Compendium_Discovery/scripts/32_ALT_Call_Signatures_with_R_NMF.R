# Run NMF on matrix based on information won from BayesNMF and subsequent analysis

# Best correlations for exposures from:
# 1) S4FB48     |   filt min100     |   K=14
# 2) P59EJM     |   filt min100     |   K=13
# 3) S4FB48     |   filt min100     |   K=14 (full TCGA)
# => all the same input matrix => SxC filt min 100
# 4) DHI7FP     |   filt min100max500   |   K=12
# 5) TOSBNB     |   filt min100 (CxS)   |   K=13
# 6) YMKE32     |   filt min100max500   |   K=13
# 7) DHI7FP     |   filt min100max500   |   K=12 (full TCGA)
#
# => Run K=13 and K=14 on SxC filt min100 and K=13 and K=12 on SxC filt min100max500

library(NMF)
library(tictoc)
library(ComplexHeatmap)

args=commandArgs(trailingOnly = TRUE)

INPUTMATRIX=args[1]
print(paste("INPUTMATRIX:", INPUTMATRIX))
K=as.numeric(args[2])
print(paste("K:", K))
CORES=as.numeric(args[3])
print(paste("CORES:", CORES))
SEED=as.numeric(args[4])
print(paste("SEED:", SEED))
NMFALG=args[5]
print(paste("NMFALG:", NMFALG))
ITER=as.numeric(args[6])
print(paste("ITER:", ITER))
OUTFOLDER=args[7]
print(paste("OUTFOLDER:", OUTFOLDER))
OUTNAME=args[8]
print(paste("OUTNAME:", OUTNAME))
dir.create(file.path(OUTFOLDER, OUTNAME), showWarnings = FALSE, recursive = TRUE)
BASENAME=file.path(OUTFOLDER, OUTNAME, paste0("5_R-NMF_K-", K, "_", OUTNAME))

## Testing purposes
# INPUTMATRIX="/Users/drews01/phd/prjcts/cnsigs2/data/3_TCGA_specific_signatures/3_CxS_Matrix_uninfprior_min100seg.rds"
# K=11
# CORES=1
# SEED=77777
# NMFALG="brunet"
# ITER=1
# OUTFOLDER="/Users/drews01/phd/prjcts/cnsigs2/data/3_TCGA_specific_signatures/5_R-NMF"
# OUTNAME="CxS_Matrix_uninfprior_min100seg"

# Load signatures
if(grepl("txt", INPUTMATRIX, ignore.case = TRUE)) {
    SxC = as.matrix( read.table(INPUTMATRIX, header = TRUE, sep = "\t", row.names = 1) )
} else {
    SxC = readRDS(INPUTMATRIX)    
}


tic()
nmfObject = NMF::nmf(SxC, rank = K, seed = SEED, nrun = ITER, method = NMFALG, .opt = paste0("p", CORES) )
toc()

saveRDS(object = nmfObject, file = paste0(BASENAME, ".rds"))

# Normalise both and save as txt and rds
if(nrow(SxC) > ncol(SxC)) {
    ## Sample by component input
    sigs = nmfObject@fit@H
    exp = nmfObject@fit@W
} else {
    ## Component by sample input
    sigs = t( nmfObject@fit@W )
    exp = t( nmfObject@fit@H )
}

rownames(sigs) = paste0("S", 1:nrow(sigs))
saveRDS(object = sigs, file = paste0(BASENAME, "_Signatures.rds"))
write.table(sigs, file = paste0(BASENAME, "_Signatures.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

sigs = t( apply(sigs, 1, function(x) x/sum(x)) )
saveRDS(object = sigs, file = paste0(BASENAME, "_Signatures_normalised.rds"))
write.table(sigs, file = paste0(BASENAME, "_Signatures_normalised.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

colnames(exp) = paste0("S", 1:ncol(exp))
exp = t( apply(exp, 1, function(x) x/sum(x)) )
saveRDS(object = exp, file = paste0(BASENAME, "_Exposures_normalised.rds"))
write.table(exp, file = paste0(BASENAME, "_Exposures_normalised.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Print Heatmaps
pdf(file = paste0(BASENAME, "_plots.pdf"), width = 12, height = 8)
Heatmap(sigs, column_title = "Normalised signatures", cluster_rows = FALSE, cluster_columns = FALSE)
Heatmap(t(exp), column_title = "Normalised exposures", cluster_rows = FALSE, cluster_columns = TRUE, show_column_names = FALSE)
dev.off()


