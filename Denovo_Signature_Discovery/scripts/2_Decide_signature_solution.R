#### Identify and extract best solution from multiple runs of BayesNMF

# You can run this script and it will simply chose the best solution according to the 
# Kullback-Leibler divergence (if option "div" chosen). Or you can manually select the most 
# representative solution based on visual inspection of network graphs of all solutions with
# optimal number of K's (rank of NMF, number of signatures). For the manual inspection, run this
# script until line 47 and then run this function to see the network graph
# with the currently chosen solution highlighted:
# testSolution(lResults$lData, lResults$allSamples, bestSample=lOptimal$dtScores$Sample[1])
# If you want this solution, simply transfer this solution to this variable and continue script:
# lOptimal$bestSample = lOptimal$dtScores$Sample[<<Whichever you find best>>]

rm(list=ls(all=TRUE))

## Libraries
library(data.table)
library(ComplexHeatmap)
library(ggplot2)
library(igraph)
library(qgraph)



## Options - TO BE EDITED BY USER
BASE="/home/drews01/SignatureDiscovery/Denovo_Signature_Discovery"
FUNCTIONS=file.path(BASE, "scripts/2_Decide_signature_solution_functions.R")
PATHTOFILES=file.path(BASE, "example_output/Out_BayesNMF_CxS_K0-43_Example_data")
OUTPUTDIR=file.path(BASE, "example_output")
# If FALSE or not numeric, then mode of K is chosen, otherwise the user-supplied K is chosen.
OVERRIDEK=10 #=> Use this for the example data
# OVERRIDEK=FALSE #=> Use this simply for the best solution from the K with the most solutions


## Load functions
source(FUNCTIONS)


#### Extract optimal solution
## Identify and load signature files from NMF
lResults = loadNMFresults(PATHTOFILES)
## Determine K
thisK = detK(lResults$lData, OVERRIDEK=OVERRIDEK)
## Identify and extract optimal solution
lOptimal = idOptimalSolution(lResults$lData, lResults$allSamples, PATHTOFILES, LOGPATTERN="log.txt", thisK = thisK, DECISION="div", OUTPUTDIR=OUTPUTDIR)

# Extract best solution for closer inspection
bestSolution = lResults$lData[[ lOptimal$bestSample ]]

###########################################
#### At this point manual inspection of solutions possible by running:
# testSolution(lResults$lData, lResults$allSamples, bestSample=lOptimal$dtScores$Sample[1])
#### If settled on a best solution, transfer the best solution to bestSolution and continue script.
# lOptimal$bestSample = <<Whatever you chose, e.g. lOptimal$dtScores$Sample[2]>>
###########################################

## Get solution for best sample (after manual selection)
bestSolution = lResults$lData[[ lOptimal$bestSample ]]

#### Plots
## Histogram of Ks
plotHist = plotHistOfKs(lResults$lData)
## How close are the signatures of the same K?
plotHeat1 = plotCosineHeat(lResults$lData, lResults$allSamples, thisK)
## Plot hairball and heatmap with best solution
plotSigs = plotHairballAndHeatmap(lResults$lData, lResults$allSamples, lOptimal$bestSample, lOptimal$thisK)
catchPlot = ViewHairball(plotSigs$hairball)

#### Combine plots into a pdf and move files from selection solution ("bestSolution") into a user-specified directory
catchTrue = makeOutputFiles(OUTPUTDIR=OUTPUTDIR, lOptimal$bestSample, lOptimal$thisK, plotHist=plotHist, plotHeat1=plotHeat1, plotSigs=plotSigs, PATHTOFILES=PATHTOFILES)
