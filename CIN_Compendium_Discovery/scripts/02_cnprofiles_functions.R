#! /usr/bin/env Rscript

library("data.table")
library("GenomicRanges")
library("ggplot2")
library("grid")
library("ggthemes")
library("Cairo")
library("cowplot")
# needed to embed Helvetica?
library("extrafont")
loadfonts(quiet = TRUE)

getChrSizes = function(CHRSIZES, chrPrefixSource = TRUE, Y = FALSE, removeChr = TRUE) {
    
    allChrSizes = fread(CHRSIZES, data.table = FALSE)
    if(isTRUE(Y)) { 
        chroms = c( seq(1:22), "X", "Y") 
    } else { 
        chroms = c( seq(1:22), "X")
    }
    
    if(isTRUE(chrPrefixSource)) { 
        relevantChrs = paste0("chr", chroms )
    } else { 
        relevantChrs = chroms
    }
    
    chrSizes = allChrSizes[ match(relevantChrs, allChrSizes$V1), ]
    if(isTRUE(removeChr)) { chrSizes$V1 = gsub("chr", "", chrSizes$V1) }
    return( chrSizes )
}

loadEvents = function( thisFileName ) {

    if( grepl( "rds", thisFileName, ignore.case = TRUE ) ) {
        theseEvents = as.data.frame( readRDS( thisFileName ) ) 
        theseEvents$nAraw = as.numeric(theseEvents$nAraw)
        theseEvents$nBraw = as.numeric(theseEvents$nBraw)
    } else {
        theseEvents = as.data.frame( fread( thisFileName ) )
    }
    eventsPerSample = split( theseEvents, theseEvents[,1] )
    return( eventsPerSample )
}

makeWGProfiles = function(thisSample, chrSizes, unrounded = TRUE, GAPS = FALSE) {
    grSample = makeGRangesFromDataFrame(thisSample, seqnames.field = "chr", start.field = "startpos", end.field = "endpos")
    
    if(isTRUE(unrounded)) {
        grSample$cn = thisSample[,7] + thisSample[,8]
        grSample$major = thisSample[,7]
        grSample$minor = thisSample[,8]
    } else {
        grSample$cn = thisSample[,5] + thisSample[,6]
        grSample$major = thisSample[,5]
        grSample$minor = thisSample[,6]
    }
    if( identical( names(seqlengths(grSample)), chrSizes$V1) ) {
        seqlengths(grSample) = chrSizes$V2
    } else {
        # X is missing and we add that.
        # To-do: if another chromosome is missing, this might not work out?
        seqlevels(grSample) = chrSizes$V1
        seqlengths(grSample) = chrSizes$V2
    }
    
    if( isTRUE(GAPS) ) {
        grGaps = gaps(grSample)
        grGaps = grGaps[strand(grGaps)=="*"]
        grGaps$cn = 0; grGaps$major = 0; grGaps$minor = 0
        grWG = sort( c(grSample, grGaps) )
        names(grWG) = seq(1:length(grWG))
        return( grWG )
    } else {
        return( grSample )
    }
}

convertForPlot = function(grWG) {
    dfWG = as.data.frame(grWG)
    
    chrSizes = seqlengths(grWG)
    toAdd = c(0, head(cumsum(as.numeric(chrSizes)), -1))
    names(toAdd) = names(seqlengths(grWG))
    
    dfWG$WGstart = dfWG$start + toAdd[ match(dfWG$seqnames, names(toAdd)) ]
    dfWG$WGend = dfWG$end + toAdd[ match(dfWG$seqnames, names(toAdd)) ]
    
    return( dfWG )
}

chrForPlot = function(grWG) { 
    chrSizes = seqlengths(grWG)
    chrNames = names(chrSizes)
    vlines = cumsum(as.numeric(chrSizes))
    textpos = as.numeric(chrSizes)/2 + c(0, head(cumsum(as.numeric(chrSizes)), -1))
    dfChrDetails = data.frame( chrNames=chrNames, vlines=vlines, textpos=textpos)
    return( dfChrDetails )
}

# Edited from: https://rpubs.com/Koundy/71792
theme_Publication <- function(base_size=14, base_family="helvetica") {
    (theme_foundation(base_size=base_size, base_family=base_family)
     + theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
             text = element_text(),
             panel.background = element_rect(colour = NA),
             plot.background = element_rect(colour = NA),
             panel.border = element_rect(colour = NA),
             axis.title = element_text(face = "bold",size = rel(1)),
             axis.title.y = element_text(angle=90,vjust =2),
             #axis.title.x = element_text(vjust = -0.2),
             axis.title.x = element_blank(),
             axis.text.y = element_text(), 
             axis.text.x = element_blank(),
             axis.line = element_line(colour="black"),
             axis.ticks.y = element_line(),
             axis.ticks.x = element_blank(),
             panel.grid.major = element_line(colour="#f0f0f0"),
             panel.grid.minor = element_blank(),
             legend.key = element_rect(colour = NA),
             # legend.position = "bottom",
             legend.position = "none",
             legend.direction = "horizontal",
             legend.key.size= unit(0.2, "cm"),
             #legend.margin = unit(0, "cm"),
             legend.spacing = unit(0, "cm"),
             legend.title = element_text(face="italic"),
             plot.margin=unit(c(10,5,5,5),"mm"),
             strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
             strip.text = element_text(face="bold")
     ))
    
}

plotCNProfile = function(dfWG, chrPlotDetails, maxVal = 5, smallPlot = TRUE, plotAlleles = FALSE, colAlt = FALSE, name ) {
    
    # Prepare arrows
    OUTLIERS=FALSE
    if( sum(dfWG$cn>maxVal) > 0 ) {
        outliers = dfWG[ dfWG$cn>maxVal,]
        outliers$lineend = c('butt')
        outliers$linejoin = c('mitre')
        OUTLIERS=TRUE
    }
    
    # main plot. add segments, either coloured or not
    if( isTRUE(colAlt) ) {
        dfWG$alt = as.factor( rep_len( c("#7570b3", "#e6ab02"), nrow(dfWG) ) )
        p = ggplot(dfWG, aes(x=WGstart, colour = dfWG$alt))+ ylim(0,5) +
            geom_segment(xend=dfWG$WGend, y=dfWG$cn, yend=dfWG$cn, size = 1.5 )
    } else {
        p = ggplot(dfWG, aes(x=WGstart))+ ylim(0,5) +
            geom_segment(xend=dfWG$WGend, y=dfWG$cn, yend=dfWG$cn, size = 1.5 )
    }
    
    # add allele-specific lines if wanted
    if( isTRUE(plotAlleles) ) {
    p = p + geom_segment(xend=dfWG$WGend, y=dfWG$major, yend=dfWG$major, colour = "red", alpha = 0.6) +
        geom_segment(xend=dfWG$WGend, y=dfWG$minor, yend=dfWG$minor, colour = "blue", alpha = 0.6)
    }
        
    # now add chr lines, change styles and names
    p = p + geom_vline(xintercept = chrPlotDetails$vlines, colour = "gray80", linetype = 2) +
        theme_Publication() + 
        ylab("Copy number") + ggtitle(name)

    # add red arrows indicating outlier segments
    if(isTRUE(OUTLIERS)) {
        p = p + geom_segment(data = outliers, aes(y = maxVal*0.95, yend = maxVal, xend = WGstart), colour = "red",
                             lineend = outliers$lineend, linejoin = outliers$linejoin, arrow = arrow(length = unit(0.1,"cm")) )
    }
    
    # If TRUE, only every second chr label will be printed
    if(isTRUE(smallPlot)) {
        p = p + geom_text(data = chrPlotDetails[c(TRUE, FALSE),], aes(y=maxVal, x=textpos, label=chrNames), 
                          size=3,vjust=-0.75, hjust=0.5, colour = "gray60")
    } else {
        p = p + geom_text(data = chrPlotDetails, aes(y=maxVal, x=textpos, label=chrNames), 
                          size=3,vjust=-0.75, hjust=0.5, colour = "gray60")
    }
    
    return( p )
}

plotPerSamples = function(eventsPerSample, BASE, OUTFOLDER, thisCancer, chrSizes, takeUnroundedVals = TRUE, 
                          plotCutoff = 5, adjustForSmall = TRUE, showAlleles = FALSE, plotAlt = FALSE,
                          device = "png", dpi = 320, print = TRUE, printGaps = FALSE) {
    
    # Sort eventsPerSample by number of events
    newOrdering = sort( sapply(eventsPerSample, function(x) nrow(x) ), decreasing = TRUE )
    eventsPerSample = eventsPerSample[ match(names(newOrdering), names(eventsPerSample)) ]
    
    # save sample-specific cn profiles
    pSamples = lapply(seq_along(eventsPerSample), function( i ) {
        
        thisSample = eventsPerSample[[i]]
        thisSample = thisSample[ ! thisSample$chr == "Y", ]
        thatName = names(eventsPerSample)[[i]]
        # thatName = as.character( thisSample[1,1] )
        
        # Get cn for whole genome
        grWG = makeWGProfiles(thisSample, chrSizes, unrounded = takeUnroundedVals, GAPS = printGaps )
        # Covnert to df and add correct coordinates
        dfWG = convertForPlot(grWG)
        # Get coordinates for chr labels and vertical lines
        chrPlotDetails = chrForPlot(grWG)
        pWG = plotCNProfile(dfWG, chrPlotDetails, name = paste( thisCancer, "-", thatName), maxVal = plotCutoff, 
                            smallPlot = adjustForSmall, plotAlleles = showAlleles, colAlt = plotAlt )
        
        if(isTRUE(print)) {
            outfile = file.path(BASE, OUTFOLDER, thisCancer, paste0( thisCancer, "_", thatName, "_CNProfile.", device) )
            if(device == "png") {
                CairoPNG(filename = outfile, dpi = dpi, width = 10, height = 6, units = "in")
            } else {
                CairoPDF(file = outfile, width = 10, height = 6)
            }
            print(pWG)
            dev.off()
        }
        
        return( pWG )
    })
    
    return( pSamples )
}

plotsPerCancer = function(pSamples, BASE = NULL, OUTFOLDER, thisCancer, plotsPerPage = 16, pageSizeMultiplier = 4 ) {
    sizeMult = plotsPerPage/pageSizeMultiplier
    pages = ceiling( length(pSamples)/plotsPerPage )
    nrow = ceiling(sqrt(plotsPerPage))
    ncol = ceiling(plotsPerPage/nrow)
    
    if( is.null(BASE) ) {
        outfile = file.path( OUTFOLDER, paste0( thisCancer, "_CNProfiles.pdf") )    
    } else {
        outfile = file.path(BASE, OUTFOLDER, paste0( thisCancer, "_CNProfiles.pdf") )
    }
    CairoPDF(file = outfile, width = sizeMult*10, height = sizeMult*6)
    for(i in 1:pages) {
        start = (i-1)*plotsPerPage+1
        end = i*plotsPerPage
        thesePlots = pSamples[start:end]
        multiPlots = plot_grid(plotlist=thesePlots, nrow = nrow, ncol = ncol)
        print(multiPlots)
    }
    dev.off()
}

plotMaster = function(BASE, INFOLDER, OUTFOLDER, CHRSIZES ) {
    chrSizes = getChrSizes(CHRSIZES)
    cancerFiles = list.files(file.path(BASE, INFOLDER))
    for (thisCancerFile in cancerFiles) {
        
        thisCancer = unlist(strsplit(thisCancerFile, "[.]"))[1]
        print(thisCancer)
        dir.create(file.path(BASE, OUTFOLDER, thisCancer), recursive = TRUE, showWarnings = FALSE)
        
        theseEvents = as.data.frame( fread(file.path(BASE, INFOLDER, thisCancerFile)) )

        eventsPerSample = split( theseEvents, theseEvents[,1] )
        
        # Plot per sample
        pSamples = plotPerSamples(eventsPerSample, BASE, OUTFOLDER, thisCancer, chrSizes)
        
        # Cancer overview pdf
        plotsPerCancer(pSamples, BASE, OUTFOLDER, thisCancer)
        
    }
    return()
}

prepareAbsoluteData = function(BASE, INFOLDER, SEGTABABS) {
    ## Read in ABSOLUTE data
    segTabAbs = fread(file.path(BASE, INFOLDER, SEGTABABS))
    
    # Take the Expected_HSCN_1 and _2
    # Reorder the data, change headers
    segTabAbsSort = segTabAbs[, c(1, 2, 3, 4, 15, 14, 11, 10)]
    names(segTabAbsSort) = paste0("V", seq(1,8))
    # change chr 23 to X
    segTabAbsSort$V2[segTabAbsSort$V2 == "23"] = "X"
    
    return( segTabAbsSort )
}


plotMasterAbsolute = function(BASE, INFOLDER, OUTFOLDER, SEGTABABS, SAMPINFO, CHRSIZES ) {
    segTabAbsSort = prepareAbsoluteData(BASE, INFOLDER, SEGTABABS)
    sampInfo = fread(file.path(BASE, INFOLDER, SAMPINFO))
    cancerTypes = unique(sampInfo$disease)
    
    chrSizes = getChrSizes(CHRSIZES)
    for (thisCancer in cancerTypes) {
        
        print(thisCancer)
        dir.create(file.path(BASE, OUTFOLDER, thisCancer), recursive = TRUE, showWarnings = FALSE)
        
        theseSamples = as.data.frame(unique(sampInfo[ sampInfo$disease == thisCancer, 1 ] ))[,1]
        shortSegs = segTabAbsSort[ segTabAbsSort$V1 %in% theseSamples, ]
        
        shortSegs$V1 = sub("^([^-]*-[^-]*-[^-]*-[^-]*).*", "\\1", shortSegs$V1)
        
        eventsPerSample = split( shortSegs, shortSegs$V1 )
        
        # Plot per sample
        pSamples = plotPerSamples(eventsPerSample, BASE, OUTFOLDER, thisCancer, chrSizes)
        
        # Cancer overview pdf
        plotsPerCancer(pSamples, BASE, OUTFOLDER, thisCancer)
        
    }
    
    return()
}


####  Special functions for multi cohort plots

loadCohorts = function(INFOLDERS, BASE, OUTFOLDER, segColName = "segVal", 
                       newColNames = c( "sample", "chr", "startpos", "endpos", "nMajor", "nMinor", "nAraw", "nBraw" ), CORES) {

    registerDoMC(CORES)
    listFolders = foreach(thisFolder=INFOLDERS, .final = function(thisFolder) setNames(thisFolder, names(INFOLDERS) ) ) %dopar% {
        
        thisFolderName = names(INFOLDERS)[ match(thisFolder, INFOLDERS) ]
        print( thisFolderName )
        
        cancerFiles = list.files(file.path(BASE, thisFolder), pattern = "raw")
        
        # Sometimes empty files with 1 bit size are present. Remove these files.
        filtFiles = sapply( cancerFiles, function(thisFile) {
            thisFileName = file.path(BASE, thisFolder, thisFile)
            if( file.size(thisFileName) < 10 ) { return( FALSE ) } else { return( TRUE ) } 
        } )
        cancerFiles = cancerFiles[filtFiles]
        
        names(cancerFiles) = do.call(rbind, strsplit(cancerFiles, "[.]"))[,1]
        listFolder = foreach(thisFile=cancerFiles, .final = function(thisFile) setNames(thisFile, names(cancerFiles) ) ) %dopar% {
            
            thisCancer = unlist(strsplit(thisFile, "[.]"))[1]
            print(thisCancer)
                    
            thisFileName = file.path(BASE, thisFolder, thisFile)
            eventsPerSample = loadEvents( thisFileName )
            dfEvents = rbindlist( eventsPerSample )

            if( segColName %in% colnames(dfEvents) ) { dfEvents$segVal = NULL }
            colnames(dfEvents) = newColNames
            dfEvents$cancer = thisCancer
            
            return( dfEvents )
        }

        dfFolder = rbindlist( listFolder )
        dfFolder$cohort = thisFolderName
        return( dfFolder )
    }

    dfFolders = as.data.frame( rbindlist( listFolders ) )
    saveRDS(listFolders, file.path( OUTFOLDER, "list_all_infolders.rds" ) )
    saveRDS(dfFolders, file.path( OUTFOLDER, "df_all_infolders.rds" ) )

    return( dfFolders )
}

# convert all possible columns to numeric and add cn and width
prepareCohortdfFolder = function( dfFolders, forLoop = c("startpos", "endpos", "nMajor", "nMinor", "nAraw", "nBraw") ) {
  
  dfFolders=as.data.frame(dfFolders)
  for(this in forLoop) {
    # avoid errors if there are factors involved
    dfFolders[, this] = as.numeric(as.character(dfFolders[ ,this]))
  }
  dfFolders=as.data.table(dfFolders)
  
  dfFolders$cn = dfFolders[, nAraw + nBraw ]
  dfFolders$width = dfFolders[, endpos - startpos ]
 
  return( dfFolders ) 
}


plotSamplesFromMultiCohorts = function(dfFolders, dfLinks, chrSizes, BASE, OUTFOLDER, CORES, printAlleles = FALSE, printAlt = FALSE) {

    registerDoMC(CORES)
    nullCatcher = foreach( i=1:nrow(dfLinks) ) %dopar% {

        thisSampleIDs = dfLinks[i,]
        thisID = thisSampleIDs[1]
        lSampleEvents = lapply(thisSampleIDs, function(ID) {
            return( dfFolders[ dfFolders$sample == ID, ] )
        })
        dfSampleEvents = as.data.frame( rbindlist( lSampleEvents ) )
        eventsPerSample = split(dfSampleEvents, dfSampleEvents$cohort)
        pSamples = plotPerSamples(eventsPerSample, BASE, OUTFOLDER, thisID, chrSizes, print = FALSE, showAlleles = printAlleles, plotAlt = printAlt ) 

        nullCatcher = plotsPerCancer(pSamples, BASE = NULL, OUTFOLDER, thisID, plotsPerPage = length(pSamples), pageSizeMultiplier = 2 )

    }

    return( NULL )
}


mergeProfiles = function(OUTFOLDER) {
    setwd(OUTFOLDER)
    system('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=merged.pdf *CNProfiles.pdf')
    return( NULL )
}

