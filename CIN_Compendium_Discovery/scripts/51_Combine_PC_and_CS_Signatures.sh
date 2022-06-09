#!/bin/bash
#!
#! SLURM job script for clust1
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#SBATCH --job-name=51_Combine_PC_and_CS_Signatures
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
echo "Pancancer sigs: $2"
echo "Pancancer exposures: $3"
echo "Path to cancer-specific files: $4"
echo "Pattern for cancer-specific sigs: $5"
echo "Pattern for cancer-specific exposures: $6"
echo "Output directory: $7"
echo "Cosine correlation threshold: $8"
echo "Cancer colours: $9"
echo "Function file: ${10}"

Rscript --vanilla ${1}/51_Combine_PC_and_CS_Signatures.R $2 $3 $4 $5 $6 $7 $8 $9 ${10}
