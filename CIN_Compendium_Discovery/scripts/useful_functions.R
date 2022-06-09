# Useful functions

l = function(x) length(x)
u = function(x) unique(x)
lu = function(x) length(unique(x))

cols2Num = function(dat, VERBOSE = TRUE) {
    
    # Two checks will be done to ensure a column can be converted into numeric
    # 1: Are there any NA's popping up after conversion? => avoid chars
    # 2: Is the column the same when converting back and forth? => avoid factors
    allCols = colnames(dat)
    for(thisCol in allCols) {
        
        x = as.character(dat[[thisCol]])
        
        test1 = sum(is.na(as.numeric(x)))
        if(test1 > 0) { next }
        
        test2 = identical(x, as.character(as.numeric(x)))
        if(! test2) { next }
        
        dat[[thisCol]] = as.numeric(x)
        if(VERBOSE) { paste("Column", thisCol, "successfully converted into numeric.") }
        
    }
    
    return(dat)   
}

doPathwayEnrichmentAnalysis = function(GENENAMES, SIGONLY = TRUE, ALPHA = 0.05) {
    require(ReactomePA)
    require(org.Hs.eg.db)
    require(clusterProfiler)
    
    # Get right identifier
    dfGenes = bitr(GENENAMES, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
    geneEntrez = unique(dfGenes$ENTREZID)
    
    ## Do the analysis, filter signif
    epa = enrichPathway(geneEntrez, organism = "human", readable = TRUE)
    
    if(SIGONLY) {
        sigEpa = data.table(epa@result[ epa@result$p.adjust < ALPHA, ])
        return(sigEpa)
    } else {
        return(epa)
    }
}


testGroupsInColumn = function(dt3, classCol, numCol, H0s, H1s) {
    
    require(data.table)
    dt3 = data.table(dt3)
    
    lTestsH0 = lapply(H0s, function(h0) { 
        
        H0 = dt3[[numCol]][ dt3[[classCol]] == h0 ]
        lTestsH1 = lapply(H1s, function(h1) { 
            
            H1 = dt3[[numCol]][ dt3[[classCol]] == h1 ]
            meanDiff = signif(mean(H1) - mean(H0), 4)
            pVal = signif(t.test(H1, H0, var.equal = FALSE)$p.value, 4)
            
            out = c(h1, h0, signif(mean(H1),4),signif(mean(H0),4), meanDiff, pVal)
            return(out)
            
        } )
        
        dtH1 = data.table(do.call(rbind, lTestsH1))
        dtH1 = cols2Num(dtH1)
        return(dtH1)
        
    } )
    
    dtH0 = rbindlist(lTestsH0)
    colnames(dtH0) = c("H1s", "H0s", "MeanH1", "MeanH0", "MeanDiff", "pVal")
    dtH0$pAdj = p.adjust(dtH0$pVal, method = "BH")
    return(dtH0)
    
}

## For graph plots (e.g. hairball plots)
#library(igraph)
#library(qgraph)

cosineSimForMatrices <- function(X, corr=FALSE){
    if(corr){ X = apply(X, 2, function(x){ x-mean(x) }) }
    denom = solve(diag(sqrt(diag(t(X)%*%X))))
    return( denom%*%(t(X)%*%X)%*%denom )
} 

plotCosineHeat = function(lData, allSamples, thisK, corThreshold = 0.9) {
    
    lFilt = list()
    for(thisSample in allSamples) {
        thisMat = lData[[ thisSample ]]
        if(nrow(thisMat) == thisK) {
            rownames(thisMat) = paste0( thisSample, "_", rownames(thisMat) )
            lFilt[[ thisSample ]] = thisMat
        }
    }
    
    
    matFilt = do.call(rbind, lFilt)
    corMat = cosineSimForMatrices(t(matFilt), corr = FALSE)
    corMat[ corMat < corThreshold ] = 0
    return( Heatmap(corMat, show_row_names = FALSE, show_column_names = FALSE, name = paste("K =", thisK)) )
    
}

makeGraphLayout = function(corMat, bestSample, thisK) {
    net = graph_from_adjacency_matrix(corMat, weighted = TRUE, mode = "undirected", diag = FALSE)
    
    # Get samples names
    theseSamples = sapply(rownames(corMat), function(thisSig) { strsplit(thisSig, "_")[[1]][1] })
    V(net)$Samples = theseSamples
    
    V(net)$color <- "#FFFFFF"
    V(net)$color[ V(net)$Samples == bestSample ] = "#228B22"
    
    # make TCGA larger
    V(net)$size = 3
    V(net)$size[ V(net)$Samples == bestSample ] = 10
    
    # label only TCGA
    V(net)$label = NA
    V(net)$label[ V(net)$Samples == bestSample ] = seq(1:thisK)
    
    # layout for spacing out vertices
    e = get.edgelist(net, names=FALSE)
    l = qgraph.layout.fruchtermanreingold(e, vcount = vcount(net),
                                          area=8*(vcount(net)^2),
                                          repulse.rad=(vcount(net)^3.1))
    
    lOut = list(net = net, l = l)
    return(lOut)
}

makeGraph = function(net, l) {
    plot(net, edge.width = 1, vertex.label.font = 2, vertex.label.color="black", layout=l)    
}



