#!/bin/bash
#!
#! SLURM job script for clust1
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#SBATCH --job-name=24_Combine_Mix_Model_Solutions
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ruben.drews@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH --output=/home/%u/data/phd/cnsigs2/data/logs/%x_slurm.%A.%a.out
#SBATCH --error=/home/%u/data/phd/cnsigs2/data/logs/%x_slurm.%A.%a.err

###############################################################
### Execute script                                         ####
###############################################################

echo -e "JobID: $SLURM_JOB_ID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"
echo ""
echo "Script path: $1"
echo "Input directory: $2"
echo "Folders for solutions: $3"
echo "Feature distribution to fit (P or G): $4"
echo "Cutoff for Gaussian mixture models: $5"
echo "Output file: $6"

Rscript --vanilla ${1}/24_Combine_Mix_Model_Solutions.R $2 $3 $4 $5 $6

## Example data:
# INPUT="/Users/drews01/phd/prjcts/cnsigs2/data/3_TCGA_specific_signatures"
# FOLDERS="2_fits_segsize,2_fits_copynumber,2_fits_changepoint,2_fits_bp10MB,2_fits_bpchrarm,2_fits_osCN"
# DISTRIBUTIONS="G,G,G,P,P,P"
# GAUSSCUTOFF=0.01
# OUTPUT="/Users/drews01/phd/prjcts/cnsigs2/data/3_TCGA_specific_signatures/2_combined_mixmodels.rds"
