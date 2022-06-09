# Calculate the sum-of-posterior matrix for the extracted copy number features

library(data.table)

args = commandArgs(trailingOnly = TRUE)
INPUTECNF=args[1]
print(paste("INPUTECNF:", INPUTECNF))
INPUTMODELS=args[2]
print(paste("INPUTMODELS:", INPUTMODELS))
OUTPUTMATRIX=args[3]
print(paste("OUTPUTMATRIX:", OUTPUTMATRIX))
UNINFPRIOR=as.logical(args[4])
print(paste("UNINFPRIOR:", UNINFPRIOR))

## For testing purposes
# INPUTECNF="/Users/drews01/data/phd/cnsigs2/data/3_TCGA_specific_signatures/1_tcga_filtered_ecnf.rds"
# INPUTMODELS="/Users/drews01/data/phd/cnsigs2/data/3_TCGA_specific_signatures/2_combined_mixmodels_merged_components.rds"
# OUTPUTMATRIX="/Users/drews01/data/phd/cnsigs2/data/3_TCGA_specific_signatures/3_SxC_Matrix_uninfprior"
# UNINFPRIOR=TRUE

allEcnf = readRDS(INPUTECNF)
allModels = readRDS(INPUTMODELS)

allFeatures = names(allModels)

lMats = lapply(allFeatures, function(thisFeature) {
  
  print(thisFeature)
  thisEcnf = allEcnf[[ thisFeature ]]
  thisModel = allModels[[ thisFeature ]]
  
  dat = as.numeric( thisEcnf[,2] )
  # We want a posterior, hence likelihood (density) times prior (weight)
  if( ncol(thisModel) == 2 ) {
    # Poisson model
    print("Poisson-based posterior")
    if(UNINFPRIOR){
        postDatUnscaled = sapply(1:nrow(thisModel), function(x) dpois(x = dat, lambda = thisModel[[x,"Mean"]]) )
    } else {
        postDatUnscaled = sapply(1:nrow(thisModel), function(x) dpois(x = dat, lambda = thisModel[[x,"Mean"]]) * thisModel[[x, "Weight"]] )    
    }
    
  } else {
    # Gaussian model
    print("Gauss-based posterior")
    if(UNINFPRIOR){
        postDatUnscaled = sapply(1:nrow(thisModel), function(x) dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]) )
    } else {
        postDatUnscaled = sapply(1:nrow(thisModel), function(x) dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]) * thisModel[[x, "Weight"]] )
    }
  }
  
  # Normalise densities to probabilities
  postDatScaled = data.frame( postDatUnscaled / rowSums(postDatUnscaled) )
  postDatScaled$Sample = thisEcnf[,1]
  matSxC = aggregate(. ~ Sample, postDatScaled, sum)
  rownames(matSxC) = matSxC$Sample
  matSxC$Sample = NULL
  matSxC = as.matrix(matSxC)
  
  # Should be sorted but just to be sure
  matSxC = matSxC[ , order(thisModel[,"Mean"]) ]
  colnames(matSxC) = paste0( thisFeature, 1:ncol(matSxC) )
  
  return(matSxC)
  
} )

allMats = do.call(cbind, lMats)

saveRDS(object = allMats, file = paste0(OUTPUTMATRIX, ".rds"))
write.table(x = allMats, file = paste0(OUTPUTMATRIX, ".txt"), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

allMats = t( allMats )
OUTPUTMATRIX = file.path(dirname(OUTPUTMATRIX), gsub("SxC", "CxS", basename(OUTPUTMATRIX)))
saveRDS(object = allMats, file = paste0(OUTPUTMATRIX, ".rds"))
write.table(x = allMats, file = paste0(OUTPUTMATRIX, ".txt"), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

