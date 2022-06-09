# Call exposures with YAPSA from a signature matrix and a sum-of-posterior matrix

args=commandArgs(trailingOnly = TRUE)

CXSMATFILE=args[1]
print(paste("CXSMATFILE:", CXSMATFILE))
SIGNATUREFILE=args[2]
print(paste("SIGNATUREFILE:", SIGNATUREFILE))
OUTFILE=args[3]
print(paste("OUTFILE:", OUTFILE))
dir.create(dirname(OUTFILE), showWarnings = FALSE, recursive = TRUE)
PLOT2FILE=as.logical(args[4])
print(paste("PLOT2FILE:", PLOT2FILE))

library(YAPSA)
if(PLOT2FILE) library(ComplexHeatmap)
if(PLOT2FILE) library(reshape2)
if(PLOT2FILE) library(data.table)
    
# ## Testing
# CXSMATFILE="/Users/drews01/data/phd/cnsigs2/data/3_TCGA_specific_signatures/3_CxS_uP_noCN.rds"
# SIGNATUREFILE="/Users/drews01/data/phd/cnsigs2/data/6_Signature_Compendium/2_Signature_Compendium_SigDefFilt.rds"
# OUTFILE="/Users/drews01/data/phd/cnsigs2/data/6_Signature_Compendium/test.pdf"
# PLOT2FILE=TRUE

V = readRDS(CXSMATFILE)
if(grepl(pattern = "rds", x = SIGNATUREFILE, ignore.case = TRUE)) {
    W = readRDS(SIGNATUREFILE)    
} else {
    W = read.table(SIGNATUREFILE, header = TRUE, sep = "\t", row.names = 1)
}


# Sanity check signature matrix
if(nrow(W) < ncol(W)) W = t(W)
# Check order of components and fix if necessary
if(! identical(rownames(W), rownames(V)) ) {
    W = W[ match(rownames(V), rownames(W)), ]
}

### YAPSA needs:
## Full matrix V        mutCatalogue        components (rows) by samples (cols)    <= HAVE
## Left matrix W        sigCatalogue        components (rows) by signature (cols)   <= HAVE
## Right matrix H       expCatalogue        signature (rows) by samples (cols)     <= WANT

# in_mutation_catalogue_df => NxM => N - Features, M - Samples => Component by Sample matrix
# in_signatures_df => NxL => N - Features, L - Signatures => Component by Signature matrix 
Hraw = as.matrix( LCD( in_mutation_catalogue_df = V, in_signatures_df = W, in_per_sample_cutoff = 0 ) )
saveRDS(object = Hraw, file = paste0(OUTFILE, "_raw.rds"))
H = t( apply(Hraw, 2, function(x) x/sum(x)) )
saveRDS(object = H, file = paste0(OUTFILE, ".rds"))

if(PLOT2FILE) {
    ## Unnormalised
    pdf(file = paste0(OUTFILE, "_raw.pdf"), width = 12 , height = 12)
    show(Heatmap(t(Hraw), show_column_names = FALSE, column_title = "YAPSA-derived Exposures"))
    
    # Boxplot of exposures per signature
    dtH = as.data.table(H)
    mH = melt(dtH)
    print(ggplot(mH, aes(x = variable, y = value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
    
    dev.off()
    
    ## Normalised
    pdf(file = paste0(OUTFILE, ".pdf"), width = 12 , height = 12)
    show(Heatmap(t(H), show_column_names = FALSE, column_title = "YAPSA-derived Exposures"))
    
    # Boxplot of exposures per signature
    dtH = as.data.table(H)
    mH = melt(dtH)
    print(ggplot(mH, aes(x = variable, y = value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
    
    dev.off()
    
}
