#!/bin/bash
#!
#! SLURM job script for clust1
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#SBATCH --job-name=53_Prepare_for_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=1:00:00
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
echo "TCGA Mixture models: $2"
echo "TCGA Signature defintions: $3"
echo "TCGA exposures on TCGA samples: $4"
echo "OV Mixture models: $5"
echo "OV Signature defintions: $6"
echo "OV exposures on TCGA samples: $7"
echo "Correlation threshold: $8"
echo "No copynumber feature present: $9"
echo "Name and path of results file: ${10}"
echo "Functions file: ${11}"

Rscript --vanilla ${1}/53_Prepare_for_analysis.R $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11}

