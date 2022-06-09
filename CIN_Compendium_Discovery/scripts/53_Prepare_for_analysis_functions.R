## Functions
extractModelsFromFlexmix = function(thisModel) {
    
    allFeatures = names(thisModel)
    DTMODEL = list()
    for(thisFeature in allFeatures) {
        
        # Get OV data from flexmix objects
        thisOld = thisModel[[ thisFeature ]]
        
        if( sum(grepl("pois", thisOld@call)) > 0 ) {
            tableOld = data.frame( parameters(thisOld) )
            colnames(tableOld) = c("Mean")
        } else {
            tableOld = data.frame( t( data.frame( parameters(thisOld) ) ) )
            colnames(tableOld) = c("Mean", "SD")
        }
        
        tableOld$Weight = prior(thisOld)
        DTMODEL[[ thisFeature ]] = tableOld[ order(tableOld$Mean), ]
        
    }
    
    return(DTMODEL)
}

getConversionMatrix = function(MODELORIGIN, MODELTARGET, plotHeatmap = FALSE, uninfPrior = TRUE, noCN = FALSE) {
    
    # Load data
    modelOrigin = readRDS(MODELORIGIN)
    modelTarget = readRDS(MODELTARGET)
    
    # Convert to list
    if("flexmix" %in% class(modelOrigin[[1]])) modelOrigin = extractModelsFromFlexmix(modelOrigin)
    if("flexmix" %in% class(modelTarget[[1]])) modelTarget = extractModelsFromFlexmix(modelTarget)
    
    # For saving posterior probabilities
    BESTOVERLAP=list()
    # Just to make sure the rownames are the same afterwards.
    ROWNAMES=list()
    allFeatures = names(modelTarget)
    if(noCN) allFeatures = allFeatures[ grep("copynumber", allFeatures, invert = TRUE) ]
    for(thisFeature in allFeatures) {
        
        # Get posterior probability for mean of OV models based on TCGA models
        originDat = modelOrigin[[ thisFeature ]]$Mean
        targetMod = modelTarget[[ thisFeature ]]
        
        if( ncol(targetMod) == 2 ) {
            # Poisson model
            originDat = round(originDat)
            if(uninfPrior) {
                # Uninformative prior. As all weights would be the same, we can just drop it as we scale later anyways.
                postDatUnscaled = sapply(1:nrow(targetMod), function(x) dpois(x = originDat, lambda = targetMod[[x,"Mean"]]) )    
            } else {
                postDatUnscaled = sapply(1:nrow(targetMod), function(x) dpois(x = originDat, lambda = targetMod[[x,"Mean"]]) *
                                             targetMod[[x,"Weight"]])
            }
            
        } else {
            # Gaussian model
            if(uninfPrior) {
                # Uninformative prior. As all weights would be the same, we can just drop it as we scale later anyways.
                postDatUnscaled = sapply(1:nrow(targetMod), function(x) dnorm(x = originDat, mean = targetMod[[x,"Mean"]], 
                                                                              sd = targetMod[[x,"SD"]]) )
            } else {
                postDatUnscaled = sapply(1:nrow(targetMod), function(x) dnorm(x = originDat, mean = targetMod[[x,"Mean"]], 
                                                                              sd = targetMod[[x,"SD"]]) * targetMod[[x,"Weight"]] )
            }
        }
        postDatScaled = data.frame( postDatUnscaled / rowSums(postDatUnscaled) )
        colnames(postDatScaled) = paste0("Target_", thisFeature, 1:ncol(postDatScaled))
        
        ROWNAMES[[thisFeature]] = paste0(thisFeature, 1:nrow(postDatScaled))
        BESTOVERLAP[[thisFeature]] = postDatScaled
    }
    
    # Join entries. Not existing entries will be filled with NAs which we set to 0.
    dtConv = rbindlist(BESTOVERLAP, fill = TRUE)
    dtConv[ is.na(dtConv) ] = 0
    mConv = as.matrix(dtConv) 
    
    rownames( mConv ) = as.vector(unlist(ROWNAMES))
    
    if(plotHeatmap) show(Heatmap( t(mConv), cluster_rows = FALSE, cluster_columns = FALSE, column_title = "Conversion Matrix"))    
    return(mConv)
    
}

liftOverSigs = function(SIGS, convMatrix) {
    
    # Load signatures
    if(grepl("txt", SIGS, ignore.case = TRUE)) {
        sigs = as.matrix( read.table(SIGS, header = TRUE, sep = "\t", row.names = 1) )
    } else {
        sigs = readRDS(SIGS)    
    }
    # We need: Signature by Component
    if(nrow(sigs) > ncol(sigs)) sigs = t(sigs)
    
    # Match order of columns
    sigs = sigs[, match(rownames(convMatrix), colnames(sigs)) ]
    
    # Convert sigs and rename components
    liftedSigs = t( sigs %*% convMatrix )
    
    return(liftedSigs)
}

plotHeatmap = function(sigCos, columnTitle, legendTitle) {
    
    colFun = colorRamp2(c(0, 1), c("white", "red"))
    return( Heatmap(sigCos, cluster_rows = FALSE, cluster_columns = FALSE, col = colFun,
                    rect_gp = gpar(col = "grey60", lwd = 1.5), column_title = columnTitle,
                    heatmap_legend_param = list( title = legendTitle, legend_height = unit(6, "cm")) ) )
    
}

compareSigs = function(liftedSigs, QUERYSIGS, convMatrix, cosineThresh = 0.85, showPlot = TRUE) {
    
    # Load signatures
    if(grepl("rds", basename(QUERYSIGS), ignore.case = TRUE)) {
        querySigs = readRDS(QUERYSIGS)
    } else {
        querySigs = read.table(QUERYSIGS, header = TRUE, sep = "\t", row.names = 1)
    }
    # We need: Component by Signature
    if(nrow(querySigs) < ncol(querySigs)) querySigs = t(querySigs)
    
    # Match order
    querySigs = querySigs[ match(gsub("Target_", "", rownames(liftedSigs)), rownames(querySigs)), ]
    
    # Do the cosine similarity comparison by vector
    sigCos = t( apply(liftedSigs, 2, function(x) cosine(x, querySigs)) )
    sigCos[ sigCos < cosineThresh ] = 0
    
    if(showPlot) {
        columnTitle = "Subject sigs (lifted; x-axis) vs query sigs (unchanged; y-axis)"
        legendTitle = "Cosine Sim"
        show( plotHeatmap(t(sigCos), columnTitle, legendTitle) )
    }
    
    return(sigCos)
    
}

compareExposures = function(BRITROCEXPOV, BRITROCEXPTCGA, cosineThresh = 0.85, showPlot = TRUE) {
    
    # Load signatures
    ovExp = readRDS(BRITROCEXPOV)
    tcgaExp = readRDS(BRITROCEXPTCGA)
    
    # Order and filt samples
    ovExp = ovExp[ rownames(ovExp) %in% rownames(tcgaExp), ]
    tcgaExp = tcgaExp[ rownames(tcgaExp) %in% rownames(ovExp), ]
    tcgaExp = tcgaExp[ match( rownames(ovExp), rownames(tcgaExp)), ]
    
    # Failsafe: If a TCGA signature has no exposure in OV samples, then add machine precision so the cosine function throws no error
    ovExp[ , colSums(ovExp) == 0 ] = ovExp[ , colSums(ovExp) == 0 ] + .Machine$double.eps
    tcgaExp[ , colSums(tcgaExp) == 0 ] = tcgaExp[ , colSums(tcgaExp) == 0 ] + .Machine$double.eps
    
    # Do the cosine similarity comparison by vector
    expCos = t( apply(ovExp, 2, function(x) cosine(x, tcgaExp)) )
    expCos[ expCos < THRESHOLD ] = 0
    
    if(showPlot) {
        columnTitle = "TCGA PC Exposures (x-axis) vs OV (y-axis)"
        legendTitle = "Cosine Sim"
        show( plotHeatmap(t(expCos), columnTitle, legendTitle) )
    }
    
    return(expCos)
    
}

prepareDataTable = function(sigCos, PREFIXCOLS = "OV", PREFIXROWS = "CS", TRESHOLD, CATEGORY) {
    
    colnames(sigCos) = paste0(PREFIXCOLS, 1:ncol(sigCos))
    rownames(sigCos) = paste0(PREFIXROWS, 1:nrow(sigCos))
    sigCos[sigCos < THRESHOLD] = 0
    sigCos[sigCos >= THRESHOLD] = 1
    dtSig = melt(round(sigCos))
    dtSig$Category = CATEGORY
    dtSig$Name = paste0(dtSig$Var1, dtSig$Var2)
    
    return(dtSig)
}

combineInformation = function(sigCos, expCos, PREFIXCOLS = "OV", PREFIXROWS = "CS", TRESHOLD = 0.85, 
                              CATSIG = "Definition", CATEXP = "Exposure", CATBOTH = "Def+Exp") {
    
    dtSig = prepareDataTable(sigCos, PREFIXCOLS = PREFIXCOLS, PREFIXROWS = PREFIXROWS, TRESHOLD, CATEGORY = CATSIG)
    dtExp = prepareDataTable(expCos, PREFIXCOLS = PREFIXCOLS, PREFIXROWS = PREFIXROWS, TRESHOLD, CATEGORY = CATEXP)
    
    dtSig$Category[ dtSig$Name %in% dtExp$Name[ dtExp$value == 1 ] & dtSig$value == 1 ] = CATBOTH
    dtExp = dtExp[ ! dtExp$Name %in% dtSig$Name[ dtSig$Category == CATBOTH ], ]
    
    dtBoth = rbind(dtSig, dtExp)
    dtBoth = dtBoth[dtBoth$value!=0,]
    
    plotOut = ggplot(dtBoth, aes(x = Var1, y = Var2, fill = Category)) + geom_tile() + scale_x_discrete(drop = FALSE) +
        scale_y_discrete(drop = FALSE) + scale_fill_brewer(palette = "Dark2") + theme_tufte(base_family = "", base_size = 18) + 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(plotOut)
}

