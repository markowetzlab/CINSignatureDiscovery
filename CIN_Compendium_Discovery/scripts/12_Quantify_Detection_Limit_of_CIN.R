
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)

ABSEXP=args[1]
print(paste("ABSEXP:", ABSEXP))
NORMEXP=args[2]
print(paste("NORMEXP:", NORMEXP))
OVSIGS=args[3]
print(paste("OVSIGS:", OVSIGS))
DETECTIONLIMIT=as.numeric(args[4])
print(paste("DETECTIONLIMIT:", DETECTIONLIMIT))
META=args[5]
print(paste("META:", META))
CLINICALDATA=args[6]
print(paste("CLINICALDATA:", CLINICALDATA))
OUTPATH=args[7]
print(paste("OUTPATH:", OUTPATH))

dir.create(OUTPATH, recursive = TRUE, showWarnings = FALSE)

# ## Intro
# BASE="/Users/drews01/phd/prjcts/cnsigs2"
# ABSEXP=file.path(BASE, "data/2_OV_signatures_on_TCGA/Export-matrix_OV_Sigs_on_TCGA-OV_rawExposures_12112019.rds")
# NORMEXP=file.path(BASE, "data/2_OV_signatures_on_TCGA/Export-matrix_OV_Sigs_on_TCGA-OV_12112019.rds")
# OVSIGS=file.path(BASE, "data/Geoffs_data/feat_sig_mat.rds")
# DETECTIONLIMIT=0.05
# META=file.path( BASE, "data/metadata/summary.ascatTCGA.penalty70.txt" )
# CLINICALDATA=file.path( BASE, "data/metadata/metadata_CNA_12K.RDS" )
# OUTPATH=file.path(BASE, "data/2_OV_signatures_on_TCGA")
# 
# DIAGNOSIS="Serous cystadenocarcinoma, NOS"
# PROJECT="TCGA-OV"

## Get samples with distinct HGSOC diagnosis
clinMetadata = readRDS(CLINICALDATA)
sampleIDs = clinMetadata$submitter_id[ clinMetadata$project_id == PROJECT & clinMetadata$primary_diagnosis == DIAGNOSIS ]

meta = fread(META)
metaFilt = meta[! is.na(meta$CNAs), ]
finalNames = metaFilt$name[ metaFilt$patient %in% sampleIDs ]


## Load exposures
absHgsoc = readRDS(ABSEXP)
absHgsoc = absHgsoc[ rownames(absHgsoc) %in% finalNames, ]

relHgsoc = readRDS(NORMEXP)
relHgsoc = relHgsoc[ rownames(relHgsoc) %in% finalNames, ]

ovSigs = readRDS(OVSIGS)

# Set everything we can detect (ie over the detection limit) to 0 and then calculate back the SxC matrix
absExpLim = absHgsoc
absExpLim[ relHgsoc > DETECTIONLIMIT ] = 0
notDetectSxC = absExpLim %*% t(YAPSA::normalize_df_per_dim(ovSigs, 2))


# Add elements up per feature
dfNotDetect = melt(notDetectSxC)
dfNotDetect$Var2 = gsub('[[:digit:]]+', '', dfNotDetect$Var2)
dtSummedNotDetect = aggregate(value ~ Var1 + Var2, dfNotDetect, sum)

dtSummedNotDetect = dtSummedNotDetect[ dtSummedNotDetect$value != 0, ]
p1 = ggplot(dtSummedNotDetect, aes( x = value)) + geom_histogram() + facet_wrap( . ~ Var2, scales = "free") + 
  theme_tufte(base_family = "ArialMT", base_size = 7) + 
  ylab("Absolute frequency") + xlab("Number of undetected CNAs")

# 95% quantile of feature distribution is our threshold
lThresh = lapply(unique(dtSummedNotDetect$Var2), function(feature) {
    THRESHOLD = as.numeric( ceiling( quantile( dtSummedNotDetect$value[ dtSummedNotDetect$Var2 == feature ], 0.95 ) ) )    
    return(paste0(feature, ": ", THRESHOLD))
} )
dtThresh = do.call(rbind, lThresh)

# Save outputs
pdf(file.path(OUTPATH, "plot_not_detectable_signal.pdf"), width = 3, height = 2); print(p1); dev.off()
ggsave(filename = file.path(OUTPATH, "plot_not_detectable_signal.svg"), p1, width = 3, height = 2)
ggsave(filename = file.path(OUTPATH, "plot_not_detectable_signal.png"), p1, width = 3, height = 2)
write.table(dtThresh, file.path(OUTPATH, "table_threshold_detectable_CIN.txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = FALSE)

