#!/bin/bash
#!
#! SLURM job script for clust1
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#SBATCH --job-name=32_ALT_Call_Signatures_with_R_NMF
#SBATCH --time=112:00:00
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
echo "Input matrix: $2"
echo "K: $3"
echo "Cores: $4"
echo "Seed: $5"
echo "NMF algorithm: $6"
echo "NMF Iterations: $7"
echo "Output folder: $8"
echo "Output basename: $9"

Rscript --vanilla ${1}/32_ALT_Call_Signatures_with_R_NMF.R $2 $3 $4 $5 $6 $7 $8 $9

