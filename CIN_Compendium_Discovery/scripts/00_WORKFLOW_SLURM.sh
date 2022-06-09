## This WORKFLOW documents the derivation CIN signatures from TCGA SNP6.0 array data.
## Started: 10/01/19
## Major overhaul: 02/12/19
## Prepared for publication: 10/02/22



#############
### Intro ###
#############

## Here only paths and scripts
BASE="/Users/drews01/SignatureDiscovery/CIN_Compendium_Discovery"
DATA="${BASE}/data"
SCRIPTS="${BASE}/scripts"


## Coming from github, please unpack the raw data
TCGA="${DATA}/rawdata/split-by-sample"
mkdir -p $TCGA
tar xzf ${DATA}/rawdata/penalty70_complete.tar.gz -C $TCGA


## Download the repository from previous work
git clone https://github.com/markowetzlab/CNsignatures ${DATA}/Macintyre2018
tar xf ${DATA}/Macintyre2018/data.tar -C ${DATA}/Macintyre2018


##############################################################
### Chapter 0: Input raw data, filter and smooth segments. ###
##############################################################

## Step 0.1: Input raw data and filter for QC
# Input: TCGA-XX-XXXX.segments.raw.txt
# Combine to one file per cancer type
# Output: one segment file per cancer type and some metadata files.

# Five filter criteria:
# a) purity > 40%
# b) purity < 100% (100% purity is most likely a mistake by ASCAT where it cannot find a solution due to low purity)
# c) oversegmentation due to low input quality (We discovered that some samples have spread out logR values; "smearing")
# d) germline mismatch (rarely mismatch between tumours and control)
# e) some input have wild logR values (mismatching BAF). These were filtered as they caused large number of CNAs in samples.
# This was discussed and confirmed in our meeting on Monday 7th of Jan 2019.
# This script filtes for a) to e) and saves to new folder.

# Please unpack "penalty70_complete.tar.gz" in the folder "data/rawdata"
INPUT="rawdata/split-by-sample"
SUMMARYFILE="metadata/Metadata_TCGA_ASCAT_penalty70.txt"
OUTPUTSEGMENTS="rawdata/split-by-cancer-type"
OUTPUTMETA="metadata"
PURITY=0.4

sbatch ${SCRIPTS}/01_merge_and_filt_raw_to_cancer_type.sh $SCRIPTS $DATA $INPUT $SUMMARYFILE $OUTPUTSEGMENTS $OUTPUTMETA $PURITY


## Step 0.2: Plot copy-number profiles [CAN BE SKIPPED]
# Input: one segmentation file per cancer
# Plots sorted by cancer type and number of breakpoints
# Output: 33 pdfs
INFOLDER="rawdata/split-by-cancer-type"
OUTFOLDER="1_cnprofiles"
CHRSIZES="metadata/hg19.chrom.sizes.txt"
FUNSCRIPT="${SCRIPTS}/02_cnprofiles_functions.R"
CORES=7
# Need to change cores in sbatch script manually

sbatch ${SCRIPTS}/02_cnprofiles.sh $SCRIPTS $DATA $INFOLDER $OUTFOLDER $CHRSIZES $FUNSCRIPT $CORES


### Step 0.3: Combine files and smooth segments
# Input: one segmentation file per cancer
# Output: one combined text and rds file
INFOLDER="${DATA}/rawdata/split-by-cancer-type"
FILTDAT="${DATA}/rawdata/combined.ascat.segments.filt"
OUTDAT="${DATA}/rawdata/combined.ascat.segments.smoothednormals"
SUMMARYFILE="${DATA}/metadata/Metadata_TCGA_ASCAT_penalty70.txt"
WIGGLE=0.1
IGNOREDELS="FALSE"
CORES=5
MEMORY="20G"

sbatch -N $CORES --mem $MEMORY ${SCRIPTS}/03_Combine_and_smooth.sh $SCRIPTS $INFOLDER $FILTDAT $OUTDAT $SUMMARYFILE $WIGGLE $IGNOREDELS $CORES




###############################################################
### Chapter 1: Determine detectable chromosomal instability ###
###############################################################

CHAPTER1FOLDER="2_OV_signatures_on_TCGA"

### Step 1.1: Calculate OV exposures
# Please make sure the github repo is downloaded and data.tar is unpacked. See code at the beginning of this script to perform these tasks automatically.
PREPATH="${DATA}/Macintyre2018/data/"
OUT="${DATA}/${CHAPTER1FOLDER}"
OVSIGS="${PREPATH}/feat_sig_mat.rds"
GEOFFPARAMS="${PREPATH}/component_parameters.rds"
OLDMETHODS="${PREPATH}/main_functions.R"
NEWMETHODS="${SCRIPTS}/main_functions.R"
META="${DATA}/metadata/Metadata_TCGA_ASCAT_penalty70.txt"
OLDCN="${DATA}/Ovarian_CN_signatures/tcga_CN_features.rds"
NEWRAW="${DATA}/rawdata/combined.ascat.segments.filt.rds"
NEWPLUS="${DATA}/rawdata/combined.ascat.segments.smoothednormals.rds"
DATE="03172022"

CORES=7
MEMORY="20G"

sbatch -N $CORES --mem $MEMORY ${SCRIPTS}/11_Compute_OV_sigs_on_TCGA.sh $SCRIPTS $PREPATH $OUT $OVSIGS $GEOFFPARAMS $CORES $OLDMETHODS $NEWMETHODS $META $OLDCN $NEWRAW $NEWPLUS $DATE


### Step 1.2: Get CIN threshold
ABSEXP="${DATA}/${CHAPTER1FOLDER}/Export-matrix_OV_Sigs_on_TCGA-OV_rawExposures_03172022.rds"
NORMEXP="${DATA}/${CHAPTER1FOLDER}/Export-matrix_OV_Sigs_on_TCGA-OV_03172022.rds"
OVSIGS="${DATA}/Macintyre2018/data/feat_sig_mat.rds"
DETECTIONLIMIT=0.05
META="${DATA}/metadata/Metadata_TCGA_ASCAT_penalty70.txt"
CLINICALDATA="${DATA}/metadata/metadata_CNA_12K.RDS"
OUTPATH="${DATA}/${CHAPTER1FOLDER}"

sbatch ${SCRIPTS}/12_Quantify_Detection_Limit_of_CIN.sh $SCRIPTS $ABSEXP $NORMEXP $OVSIGS $DETECTIONLIMIT $META $CLINICALDATA $OUTPATH


### Step 1.3: Detect CIN
META="${DATA}/metadata/Metadata_TCGA_ASCAT_penalty70.txt"
# Samples with CNAs > 19 (Minimum of 20 CNAs) have detectable CIN
CINTHRESHOLD=19
NEWPLUS="${DATA}/rawdata/combined.ascat.segments.smoothednormals.rds"
OUTFOLDER="${DATA}/3_Pancancer_Signatures"

sbatch ${SCRIPTS}/13_Detect_CIN.sh $SCRIPTS $META $CINTHRESHOLD $NEWPLUS $OUTFOLDER




#####################################
### Chapter 2: Fit mixture models ###
#####################################


CHAPTER2FOLDER="3_Pancancer_Signatures"

### Step 2.1: Extract copy-number features
# Prepare: Mixture models
SOURCE="${SCRIPTS}/main_functions.R"
INPUTFILE="${DATA}/${CHAPTER2FOLDER}/0_TCGA_Segments_dCIN.rds"
OUTPATH="${DATA}/${CHAPTER2FOLDER}"
CORES=7
PREPATH="${DATA}/Macintyre2018/data/"
RMNORM=TRUE

sbatch -N $CORES ${SCRIPTS}/21_ECNF_TCGA.sh $SCRIPTS $SOURCE $INPUTFILE $OUTPATH $CORES $PREPATH $RMNORM



### Step 2.2: Plot feature distributions
# Control: To check what to expect from the mixture models
SCRIPTMARKDOWN="${SCRIPTS}/22_Plot_ECNF.Rmd"
INPUTECNF="${DATA}/${CHAPTER2FOLDER}/1_tcga_filtered_ecnf.rds"
OUTPATH="${DATA}/${CHAPTER2FOLDER}"

sbatch ${SCRIPTS}/22_Plot_ECNF.sh $SCRIPTMARKDOWN $INPUTECNF $OUTPATH



### Step 2.3A: Fit mixture models (Finite Mixture Models in R)
###
### Important notes on the different mixture model fits!
###
### SEGSIZE / CHANGEPOINT / COPYNUMBER (Gaussian mixture models)
### The three features segsize, changepoint and copynumber have been modelled with a Dirichlet-process Gaussian mixture model using variational inference.
### The results are in scripts/viMixtures and need to be copied into ${CHAPTER2FOLDER}/2_fits_10k_segsize/ and ${CHAPTER2FOLDER}/2_fits_5k_changepoint/.
### The result file is a csv file containing all mixtures. The summary script "24_Combine_Mix_Model_Solutions" is filtering the mixtures
### by a minimum weight criterion and these components are then being used in script XX to call the sum-of-posterior(-probabilities) for
### each event summed over for each sample. The lines underneath are running mclust for segsize, changepoint and copynumber to illustrate
### how they are not able to fit the data reasonably well. Hence this code is uncommented but present.
###
### BPCHRARM / BP10MB / OSCN (Poisson mixture models)
### The features bpchrarm, bp10MB and osCN are still modelled by Poisson mixtures from the flexmix packages. Here, I also deviated from
### Geoff's original plan by running the fits 100 times, let K decide by the mode of all 100 runs and then take the run with the best BIC
### value. This is done in script "24_Combine_Mix_Model_Solutions.sh".


### Intro
# One job per feature distribution
PATHMAIN="${SCRIPTS}/main_functions.R"
INPUTECNF="${DATA}/${CHAPTER2FOLDER}/1_tcga_filtered_ecnf.rds"
OUTPATH="${DATA}/${CHAPTER2FOLDER}"

### Uncommented but present. See explanation above.
# ## segsize, changepoint and copynumber (Normal distribution)
# ## Showing the inability to find a consistent number of clusters with finite Gaussian mixture models.
# FEATURES="segsize changepoint copynumber"
#
# for FEATURE in $FEATURES; do
#
#     # Caveat: The range of scanned components should not be larger than the iterations. Otherwise overlap of seeds and thus non-unique solutions.
#     MINCOMP=1
#     MAXCOMP=30
#     ITER=100
#     CORES=10
#     MEMORY="30G"
#
#     sbatch -N $CORES --mem $MEMORY ${SCRIPTS}/23_Fit_MixtureModels_Normals.sh $SCRIPTS $INPUTECNF $OUTPATH $FEATURE $MINCOMP $MAXCOMP $ITER $CORES
#
# done

## bp10MB (Poisson distribution)
# Previous successful prior: 0.001 and with K=3
FEATURE="bp10MB"
FITDIST="pois"
MINCOMP=2
MAXCOMP=10
PRIOR=0.001
NITER=1000
ITERS=100

DECISION="BIC"
for ITER in $(seq 1 $ITERS); do

    UUID=$(cat /dev/urandom | tr -dc 'A-Z0-9' | fold -w 6 | head -n 1)
    OUTPATH="${DATA}/${CHAPTER2FOLDER}/2_fits_bp10MB/${UUID}"
    mkdir -p $OUTPATH
    sbatch ${SCRIPTS}/23_Fit_MixtureModels.sh $SCRIPTS $PATHMAIN $INPUTECNF $OUTPATH $FEATURE $FITDIST $MINCOMP $MAXCOMP $PRIOR $NITER $DECISION $ITER

done

## osCN (Poisson distribution)
# Previous successful prior: 0.001 and with K=4
FEATURE="osCN"
FITDIST="pois"
MINCOMP=2
MAXCOMP=10
PRIOR=0.001
NITER=1000
ITERS=100

DECISION="BIC"
for ITER in $(seq 1 $ITERS); do

    UUID=$(cat /dev/urandom | tr -dc 'A-Z0-9' | fold -w 6 | head -n 1)
    OUTPATH="${DATA}/${CHAPTER2FOLDER}/2_fits_osCN/${UUID}"
    mkdir -p $OUTPATH
    sbatch ${SCRIPTS}/23_Fit_MixtureModels.sh $SCRIPTS $PATHMAIN $INPUTECNF $OUTPATH $FEATURE $FITDIST $MINCOMP $MAXCOMP $PRIOR $NITER $DECISION $ITER

done

## bpchrarm (Poisson distribution)
# Previous successful prior: 0.001 and with K=6
FEATURE="bpchrarm"
FITDIST="pois"
MINCOMP=2
MAXCOMP=10
PRIOR=0.001
NITER=3000
ITERS=100

DECISION="BIC"
for ITER in $(seq 1 $ITERS); do

    UUID=$(cat /dev/urandom | tr -dc 'A-Z0-9' | fold -w 6 | head -n 1)
    OUTPATH="${DATA}/${CHAPTER2FOLDER}/2_fits_bpchrarm/${UUID}"
    mkdir -p $OUTPATH
    sbatch ${SCRIPTS}/23_Fit_MixtureModels.sh $SCRIPTS $PATHMAIN $INPUTECNF $OUTPATH $FEATURE $FITDIST $MINCOMP $MAXCOMP $PRIOR $NITER $DECISION $ITER

done


### Step 2.3B: Fit mixture models (Variational Inference DPGMM in Python)
## Extract feature distributions
INPUTECNF="${DATA}/${CHAPTER2FOLDER}/1_tcga_filtered_ecnf.rds"
OUTPATH="${DATA}/${CHAPTER2FOLDER}"

Rscript --vanilla ${SCRIPTS}/23B_Preparation_Extract_Distributions.R $INPUTECNF $OUTPATH


## Run viMixtures (currently not opmtised for HPC usage)
# Use the spec-file.txt in scripts/viMixtures/ to create the proper
# conda environment to run the Python script
# conda create --name vidpgmm --file ${SCRIPTS}/viMixtures/spec-file.txt
# Then afterwards, run:
# conda activate vidpgmm
INPUTSS="${DATA}/${CHAPTER2FOLDER}/1_tcga_segmentsize_dist.csv"
OUTPUT="${DATA}/${CHAPTER2FOLDER}/2_VI_DPGMMs"
mkdir -p $OUTPUT
python ${SCRIPTS}/viMixtures/inferMixtures.py $INPUTSS > ${OUTPUT}/tcga_changepoint_4_5000.csv

INPUTCP="${DATA}/${CHAPTER2FOLDER}/1_tcga_changepoint_dist.csv"
python ${SCRIPTS}/viMixtures/inferMixtures.py $INPUTCP > ${OUTPUT}/tcga_segsize_4_10000.csv


### Step 2.4: Chose optimal mixture models
# Input: 5 directories
# Intermediate: For Poisson runs, chose best K by taking mode of all Ks and then decide for best fit amongst this K
# Output: one object with six models (only data table with means+weights or means+sd+weights)
INPUT="${DATA}/${CHAPTER2FOLDER}"
FOLDERS="2_fits_10k_segsize,2_fits_5k_changepoint,2_fits_bp10MB,2_fits_bpchrarm,2_fits_osCN"
DISTRIBUTIONS="G,G,P,P,P"
GAUSSCUTOFF=0.01
OUTPUT="${DATA}/${CHAPTER2FOLDER}/2_combined_mixmodels_5k.rds"

sbatch ${SCRIPTS}/24_Combine_Mix_Model_Solutions.sh $SCRIPTS $INPUT $FOLDERS $DISTRIBUTIONS $GAUSSCUTOFF $OUTPUT



### Step 2.5: Merge compontents
# Now manuallly inspect mixture model as we have found that strongly overlapping components introduce artefacts in the
# signatures, ie one component in one sig and the other component in another sig although they are almost identical.
# Merging done by weighted average of mean and sum of weights. For sd I have drawn 100k points from the original distributions and
# then estimated sd from there.
# Merging criteria:
# - Closer than 1MB for medium- and large-size elements OR roughly 10% of size OR strong overlap due to large SDs.
# - For CN roughly closer than 0.1.
# See Supplementary Table 15 and then manually replace components from the
# original mixture components with the newly identified merged mixture compontents

# Run this script and the mixture components get automatically replaced to the version used in the paper.
INPUT="${DATA}/${CHAPTER2FOLDER}/2_combined_mixmodels_5k.rds"
OUTPUT="${DATA}/${CHAPTER2FOLDER}/2_combined_mixmodels_merged_components.rds"
Rscript --vanilla ${SCRIPTS}/helper_scripts/Replace_Mixture_Components.R $INPUT $OUTPUT




#######################################################
### Chapter 3: Derive pan-cancer signatures of CNAs ###
#######################################################


CHAPTER2FOLDER="3_Pancancer_Signatures"

### Step 3.1: Calculate sum-of-posterior matrix
INPUTECNF="${DATA}/${CHAPTER2FOLDER}/1_tcga_filtered_ecnf.rds"
INPUTMODELS="${DATA}/${CHAPTER2FOLDER}/2_combined_mixmodels_merged_components.rds"
OUTPUTMATRIX="${DATA}/${CHAPTER2FOLDER}/3_SxC_uninfPrior"
UNINFPRIOR="TRUE"

sbatch ${SCRIPTS}/31_Calculate_SxC_Matrix.sh $SCRIPTS $INPUTECNF $INPUTMODELS $OUTPUTMATRIX $UNINFPRIOR



### Step 3.2: Call pan-cancer signatures
# This part is done offline on a computer with a GPU by using Python and Pytorch.
# See the "Denovo_Signature_Discovery" part of the repo for more details.
NMFPATH="/Users/drews01/SignatureDiscovery/Denovo_Signature_Discovery"
SCRIPTPATH="${NMFPATH}/SignatureAnalyzer-GPU/SignatureAnalyzer-GPU.py"
INPUTMATRIX="${DATA}/${CHAPTER2FOLDER}/3_CxS_uninfPrior.txt"
MAXITER=100000
K0=50
OUTPUTDIR="${DATA}/${CHAPTER2FOLDER}/4_NMF_CxS_K0-${K0}"
REGULARISATIONS="--prior_on_W L1 --prior_on_H L1"

## If you want to seed a different H (signature) matrix.
# K0=17
# OUTPUTDIR="/Users/drews01/phd/prjcts/cnsigs2/data/${CHAPTER2FOLDER}/xx_sigs_from_joint_pccs_seed_mat-K0-${K0}"
# PRIORH="/Users/drews01/phd/prjcts/cnsigs2/data/5_Cancer_Signature_Compendium/2_Signature_Compendium_onlySigs.txt"
# PRIORW="/Users/drews01/phd/prjcts/cnsigs2/data/5_Cancer_Signature_Compendium/2_Exposure_Compendium_onlySigs.txt"
# REGULARISATIONS="--prior_on_W L1 --prior_on_H L1 --prior_mat_H ${PRIORH} --prior_mat_W ${PRIORW}"

ITERS=12
# Every independent run takes up roughl 545MB graphics memory.
# ~9400 free ram (not OS) => 17 parallel runs possible
PARALLEL=16

# Thanks to the python script I have to go the directory from where the outputdir is then added.
mkdir -p $OUTPUTDIR
cd $OUTPUTDIR;
conda activate bayesnmf
#export PYTHONPATH="/home/drews01/anaconda3/envs/bayesnmf/lib/python3.6/site-packages:/home/drews01/anaconda3/envs/bayesnmf/lib/python3.7/site-packages:$PYTHONPATH"
for ITER in $(seq 1 $ITERS); do

    echo $ITER

    for PAR in $(seq 1 $PARALLEL); do

        UUID=$(cat /dev/urandom | tr -dc 'A-Z0-9' | fold -w 6 | head -n 1)
        OUTPATH="${UUID}"
        mkdir -p $OUTPATH
        python $SCRIPTPATH --data $INPUTMATRIX --max_iter=$MAXITER --output_dir $OUTPATH $REGULARISATIONS --labeled --K0 $K0 > ${OUTPATH}/${UUID}_log.txt &

    done

    sleep 4
    PYRUN=`nvidia-smi | grep -c python`
    while [[ $PYRUN != 0 ]]; do
        sleep 2
        PYRUN=`nvidia-smi | grep -c python`
    done

done



### ALTERNATIVE Step 3.2: Call pan-cancer signatures with R NMF
# Note: As CPU-based NMF is insanely slow on large matrices, you need to know the K beforehand! BayesNMF is again a very good indicator for that.
# Input: SxC matrix
# Output: RDS object with signature and exposure matrix for a give K
CORES=15
SEED=77777
NMFALG="brunet"
ITER=15
OUTFOLDER="${DATA}/${CHAPTER2FOLDER}/5_R-NMF"

INPUTMATRIX="${DATA}/${CHAPTER2FOLDER}/3_SxC_uP_noCN.txt"
MEMORY="30G"
NAME="SxC_uP_noCN"
K=9
sbatch -N $CORES --mem $MEMORY ${SCRIPTS}/32_ALT_Call_Signatures_with_R_NMF.sh $SCRIPTS $INPUTMATRIX $K $CORES $SEED $NMFALG $ITER $OUTFOLDER $NAME
K=10
sbatch -N $CORES --mem $MEMORY ${SCRIPTS}/32_ALT_Call_Signatures_with_R_NMF.sh $SCRIPTS $INPUTMATRIX $K $CORES $SEED $NMFALG $ITER $OUTFOLDER $NAME

INPUTMATRIX="${DATA}/${CHAPTER2FOLDER}/3_CxS_uP_noCN.txt"
MEMORY="90G"
NAME="CxS_uP_noCN"
K=9
sbatch -N $CORES --mem $MEMORY ${SCRIPTS}/32_ALT_Call_Signatures_with_R_NMF.sh $SCRIPTS $INPUTMATRIX $K $CORES $SEED $NMFALG $ITER $OUTFOLDER $NAME
K=10
sbatch -N $CORES --mem $MEMORY ${SCRIPTS}/32_ALT_Call_Signatures_with_R_NMF.sh $SCRIPTS $INPUTMATRIX $K $CORES $SEED $NMFALG $ITER $OUTFOLDER $NAME



### Step 3.3: Decide for optimal pan-cancer signature
# Input: Folders from BayesNMF
# Output: RDS and txt objects of best solution
PATHTOFILES="${DATA}/${CHAPTER2FOLDER}/4_NMF_CxS_K0-50"
SIGMATPATTERN="W.txt"
EXPMATPATTERN="H.txt"
LOGPATTERN="log.txt"
corThreshold=0.85
#OUTPUTDIR="${DATA}/3_TCGA_specific_signatures"
OUTPUTDIR=$PATHTOFILES
TRANSPOSED=FALSE
DECISION="div"
WHICHK="NULL"
METADATA="${DATA}/metadata/Metadata_TCGA_ASCAT_penalty70.txt"
CANCERCOLS="${DATA}/metadata/TCGA_colour_scheme_Lydia.txt"

sbatch ${SCRIPTS}/33_Decide_Pancancer_Signature.sh $SCRIPTS $PATHTOFILES $SIGMATPATTERN $EXPMATPATTERN $LOGPATTERN $corThreshold $OUTPUTDIR $TRANSPOSED $DECISION $WHICHK $METADATA $CANCERCOLS


########################################################################################################################
#### MANUAL: DECIDE ON WHICH SET OF PANCANCER SIGNATURES TO FOCUS ON!                                               ####
#### What does it mean? Have a look at page 5 of the summary pdf (the graph with the dots). We chose a solution     ####
#### which had a signature from each cluster. It might not be the mathematically best solution but it represents    ####
#### the best out of all runs with K=9.                                                                             ####
#### For more accessible code and better documentation, see folder "Denovo_Signature_Discovery".                    ####
########################################################################################################################




############################################################
### Chapter 4: Derive cancer-specific signatures of CNAs ###
############################################################

CHAPTER2FOLDER="3_Pancancer_Signatures"
CHAPTER4FOLDER="4_Cancer_specific_Signatures"

### Step 4.1: Decide which cancers are suitable and split SxC matrix
METADATA="${DATA}/metadata/Metadata_TCGA_ASCAT_penalty70.txt"
MINSAMPS=100
IGNOREBRCASPLIT=TRUE
SXCMATRIX="${DATA}/${CHAPTER2FOLDER}/3_SxC_uninfPrior.rds"
OUTPUTDIR="${DATA}/${CHAPTER4FOLDER}"
mkdir -p $OUTPUTDIR

sbatch ${SCRIPTS}/41_Find_Suitable_Cancers_and_Split_SxC.sh $SCRIPTS $METADATA $MINSAMPS $IGNOREBRCASPLIT $SXCMATRIX $OUTPUTDIR


### Step 4.2: Call cancer-specific signatures
# This part is done offline on a computer with a GPU by using Python and Pytorch.
# Input: Directory with folders containing SxC matrices
# Output: x runs of BayesNMF
NMFPATH="/Users/drews01/SignatureDiscovery/Denovo_Signature_Discovery"
SCRIPTPATH="${NMFPATH}/SignatureAnalyzer-GPU/SignatureAnalyzer-GPU.py"
PATHTOFILES="${DATA}/${CHAPTER4FOLDER}"
# BayesNMF and ARD suggest V=WH => W being the dictionary and H the activation matrix
MATRIX="3_CxS"
ALLCANCERS=`ls $PATHTOFILES | grep -Ev ".txt|.rds"`

for thisCancer in $ALLCANCERS; do

    echo $thisCancer

    INPUTMATRIX=`ls ${PATHTOFILES}/${thisCancer}/${MATRIX}_${thisCancer}.txt`
    MAXITER=100000
    OUTPUTDIR="${PATHTOFILES}/${thisCancer}/BayesNMF_runs"
    mkdir -p $OUTPUTDIR
    REGULARISATIONS="--prior_on_W L1 --prior_on_H L1"

    ITERS=12
    # Every independent run takes up roughly 450MB graphics memory.
    PARALLEL=17

    # Thanks to the python script I have to go the directory from where the outputdir is then added.
    cd $OUTPUTDIR
    conda activate bayesnmf

    for ITER in $(seq 1 $ITERS); do

        echo $ITER

        for PAR in $(seq 1 $PARALLEL); do

            UUID=$(cat /dev/urandom | tr -dc 'A-Z0-9' | fold -w 6 | head -n 1)
            OUTPATH="${UUID}"
            mkdir -p $OUTPATH
            python $SCRIPTPATH --data $INPUTMATRIX --max_iter=$MAXITER --output_dir $OUTPATH $REGULARISATIONS --labeled > ${OUTPATH}/${UUID}_log.txt &

        done

        sleep 4
        PYRUN=`nvidia-smi | grep -c python`
        while [[ $PYRUN != 0 ]]; do
            sleep 2
            PYRUN=`nvidia-smi | grep -c python`
        done

    done

done


### Step 4.3: Decide on cancer-specific signatures
PATHTOFILES="${DATA}/${CHAPTER4FOLDER}"
ALLCANCERS=`ls $PATHTOFILES | grep -Ev ".txt|.rds"`
SIGMATPATTERN="W.txt"
EXPMATPATTERN="H.txt"
LOGPATTERN="log.txt"
corThreshold=0.85
TRANSPOSED=FALSE
DECISION="div"
ALLK="1 2"
METADATA="${DATA}/metadata/Metadata_TCGA_ASCAT_penalty70.txt"
CANCERCOLS="${DATA}/metadata/TCGA_colour_scheme.txt"

for WHICHK in $ALLK; do

    for thisCancer in $ALLCANCERS; do

        echo $thisCancer

        PATHTOFILES="${DATA}/${CHAPTER4FOLDER}/${thisCancer}/BayesNMF_runs"
        OUTPUTDIR="${DATA}/${CHAPTER4FOLDER}/${thisCancer}"

        sbatch ${SCRIPTS}/33_Decide_Pancancer_Signature.sh $SCRIPTS $PATHTOFILES $SIGMATPATTERN $EXPMATPATTERN $LOGPATTERN $corThreshold $OUTPUTDIR $TRANSPOSED $DECISION $WHICHK $METADATA $CANCERCOLS

    done

done

########################################################################################################################
#### MANUAL: DECIDE ON WHICH SET OF CANCER-SPECIFIC SIGNATURES TO FOCUS ON!                                         ####
#### What does it mean? Have a look at page 5 of the summary pdf (the graph with the dots). We chose a solution     ####
#### which had a signature from each cluster. It might not be the mathematically best solution but it represents    ####
#### the best out of all runs with that K. Sometimes chosing the second best K to reach more stable solutions.      ####
#### Have a look at the "all_Solutions" file which documents the order of best solutions.                           ####
#### For more accessible code and better documentation, see folder "Denovo_Signature_Discovery".                    ####
########################################################################################################################




###################################################################################################
### Chapter 5: Create signature compendium by merging pan-cancer and cancer-specific signatures ###
###################################################################################################

CHAPTER1FOLDER="2_OV_signatures_on_TCGA"
CHAPTER2FOLDER="3_Pancancer_Signatures"
CHAPTER4FOLDER="4_Cancer_specific_Signatures"
CHAPTER5FOLDER="5_Signature_Compendium"

## Step 5.1: Combine pan-cancer and cancer-specific signatures
# Input:
# Pancancer signatures
# Cancer-specific folders with BayesNMF output files in it
# Output:
# Plots comparing pancancer and cancer-specific signatures
# A set of cancer-related signatures
ID="KRJ5F9" ## This needs replacement if you run your own signatures.
PANCANCERSIGS="${DATA}/${CHAPTER2FOLDER}/4_Signatures_${ID}_normalised.rds"
PANCANCEREXP="${DATA}/${CHAPTER2FOLDER}/4_Exposures_${ID}.rds"
PATHCANCERSPECIFIC="${DATA}/${CHAPTER4FOLDER}"
PATTERNSIGS="4_Signatures"
PATTERNEXP="4_Exposures"
OUTDIR="${DATA}/${CHAPTER5FOLDER}"
# To determine a sensible threshold, please have a look at Cosine_Thresh_CS_PC_Merge_Analysis.R which simulated signature compendiums and
# then estimates a threshold.
COSINETHRESH=0.74
CANCERCOLS="${DATA}/metadata/TCGA_colour_scheme.txt"
SOURCEFUNCTIONS="${SCRIPTS}/51_Combine_PC_and_CS_Signatures_functions.R"

sbatch ${SCRIPTS}/51_Combine_PC_and_CS_Signatures.sh $SCRIPTS $PANCANCERSIGS $PANCANCEREXP $PATHCANCERSPECIFIC $PATTERNSIGS \
    $PATTERNEXP $OUTDIR $COSINETHRESH $CANCERCOLS $SOURCEFUNCTIONS



## Step 5.2: Call exposures on TCGA samples for CSS and Compendium signatures
CXSMATFILE="${DATA}/${CHAPTER2FOLDER}/3_CxS_uninfPrior.rds"
PLOT2FILE=TRUE

# Whole Signature Compendium
SIGNATUREFILE="${DATA}/${CHAPTER5FOLDER}/2_Signature_Compendium_Cosine-${COSINETHRESH}.rds"
OUTFILE="${DATA}/${CHAPTER5FOLDER}/3_Exposures_Signature_Compendium_Cosine-${COSINETHRESH}"
sbatch ${SCRIPTS}/52_Call_Exposures_from_Signatures.sh $SCRIPTS $CXSMATFILE $SIGNATUREFILE $OUTFILE $PLOT2FILE

# Only CSS
# SIGNATUREFILE="${DATA}/${CHAPTER5FOLDER}/2_Signatures_CSS_Cosine-0.85.rds"
# SIGNATUREFILE="${DATA}/${CHAPTER4FOLDER}/OV"
# OUTFILE="${DATA}/${CHAPTER2FOLDER}/OV/5_Exposures_OV"
# sbatch ${SCRIPTS}/52_Call_Exposures_from_Signatures.sh $SCRIPTS $CXSMATFILE $SIGNATUREFILE $OUTFILE $PLOT2FILE



## Step 5.3: Compare compendium with OV signatures for first validation
# This script compares the compendium signatures to the OV signatures based on their signature definitions and their exposure correlations
# on TCGA-OV samples. The TCGA-OV cohort was preferred over the BRITROC cohort as there is SNP array data, larger cohort and have been
# analysed with the same CN feature extraction methods.
# It also renames the signatures based on the number of exposed samples and saves the output in the results folder.
MIXMODELSTCGA="${DATA}/${CHAPTER2FOLDER}/2_combined_mixmodels_merged_components.rds"
SIGDEFSTCGA="${DATA}/${CHAPTER5FOLDER}/2_Signature_Compendium_Cosine-${COSINETHRESH}.rds"
EXPTCGATCGA="${DATA}/${CHAPTER5FOLDER}/3_Exposures_Signature_Compendium_Cosine-${COSINETHRESH}.rds"

MIXMODELSOV="${DATA}/Macintyre2018/data/component_parameters.rds"
SIGDEFSOV="${DATA}/Macintyre2018/data/feat_sig_mat.rds"
EXPOVTCGA="${DATA}/${CHAPTER1FOLDER}/Export-matrix_OV_Sigs_on_TCGA-OV_03172022.rds"

THRESHOLD=0.85
NOCN=TRUE
RESULTSFILE="${BASE}/results/Signature_Compendium_v5_Cosine-${COSINETHRESH}"
FUNCTIONS="${SCRIPTS}/53_Prepare_for_analysis_functions.R"

sbatch ${SCRIPTS}/53_Prepare_for_analysis.sh $SCRIPTS $MIXMODELSTCGA $SIGDEFSTCGA $EXPTCGATCGA $MIXMODELSOV $SIGDEFSOV $EXPOVTCGA $THRESHOLD $NOCN $RESULTSFILE $FUNCTIONS
