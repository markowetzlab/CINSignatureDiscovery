# Combine the mix model solutions to one R object

# Assumptions:
## If distribution "G", then a csv or other text file in this folder with four columns Index, Means, Variance, Weight.
## If distribution "P", then many subfolders with flexmix runs. Chose the mode of K from multiple runs via the *simple.txt files
## and chose the best solution which has the best BIC from all solutions with K=Koptimal.

library(data.table)
library(flexmix)

args = commandArgs(trailingOnly = TRUE)
INPUT=args[1]
print(paste("INPUT:", INPUT))
FOLDERS=args[2]
print(paste("FOLDERS:", FOLDERS))
DISTRIBUTIONS=args[3]
print(paste("DISTRIBUTIONS:", DISTRIBUTIONS))
GAUSSCUTOFF=as.numeric( args[4] )
print(paste("GAUSSCUTOFF:", GAUSSCUTOFF))
OUTPUT=args[5]
print(paste("OUTPUT:", OUTPUT))

## Testing
# INPUT="/Users/drews01/phd/prjcts/cnsigs2/data/3_TCGA_specific_signatures"
# FOLDERS="2_fits_segsize,2_fits_copynumber,2_fits_changepoint,2_fits_bp10MB,2_fits_bpchrarm,2_fits_osCN"
# DISTRIBUTIONS="G,G,G,P,P,P"
# GAUSSCUTOFF=0.01
# OUTPUT="/Users/drews01/phd/prjcts/cnsigs2/data/3_TCGA_specific_signatures/2_combined_mixmodels.rds"

getmode = function(allKs) {
    uniqKs = unique(allKs)
    uniqKs[ which.max( tabulate( match( allKs, uniqKs ) ) ) ]
}


allFolders = strsplit(FOLDERS, ",")[[1]]
allDists = strsplit(DISTRIBUTIONS, ",")[[1]]

thisOut = list()

for(thisIter in 1:length(allFolders)) {

    thisFolder = allFolders[thisIter]
    splitName = strsplit( thisFolder, "_" )[[1]]
    thisFeature = splitName[ length( splitName ) ]
    thisDist = allDists[thisIter]
    print(thisFolder)
    print(thisFeature)
    print(thisDist)


    if( thisDist == "G" ) {

        thisFile = list.files(path = file.path(INPUT, thisFolder), pattern = ".csv", full.names = TRUE)
        if(length(thisFile) > 1) { warning(paste("This feature has multiple csv files in its folder. Chose the first one:", thisFeature)) }
        theseParams = fread(thisFile)
        filtParams = theseParams[ theseParams$Weight > GAUSSCUTOFF ]
        filtParams$SD = sqrt(filtParams$Variance)
        filtParams$Index = NULL
        filtParams$Variance = NULL
        thisOut[[ thisFeature ]] = filtParams[ ,c("Mean", "SD", "Weight") ]

    } else if( thisDist == "P" ) {

        # Step 1: Decide on K by mode
        theseFiles = list.files(path = file.path(INPUT, thisFolder), pattern = ".txt", recursive = TRUE, full.names = TRUE)
        allKs = as.numeric( sapply(theseFiles, function(thisFile) { nrow(fread(thisFile)) }) )
        thisK = getmode(allKs)
        print(paste("Mode is", thisK))

        # Step 2: Take the model with the best fit for this given K
        # Important: Every run must be finished, so every run has to have a rds and a txt file. Otherwise this part will be faulty.
        thisKFiles = list.files(path = file.path(INPUT, thisFolder), pattern = "rds", ignore.case = TRUE,
                                recursive = TRUE, full.names = TRUE)[ allKs == thisK ]
        bestSolutionBIC = Inf
        bestSolution = NA
        for(thisFile in thisKFiles) {

            thisModel = readRDS(thisFile)
            thisBIC = abs(BIC(thisModel))
            if(thisBIC < bestSolutionBIC) {
                print(paste("Currently this model is the best fitting model:", thisFile))
                bestSolutionBIC = thisBIC
                bestSolution = thisModel
            }

        }

        dtSolution = data.table(Mean = parameters(bestSolution), Weight = prior(bestSolution))
        dtSolution = dtSolution[ order(dtSolution$Mean) ]
        thisOut[[ thisFeature ]] = dtSolution

    }

}

str(thisOut)
saveRDS(object = thisOut, file = OUTPUT)
