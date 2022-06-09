#!/bin/bash
#!
#! SLURM job script for clust1
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#SBATCH --job-name=03_Combine_and_smooth
#SBATCH --time=1:00:00
#SBATCH --mail-type=END,FAIL
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
echo "Input folder: $2"
echo "Combined output file: $3"
echo "Combined and smoothed output file: $4"
echo "Kerstin's summary file: $5"
echo "Wiggle room for smoothing: $6"
echo "Ignore deletions when smoothing: $7"
echo "Cores: $8"

Rscript --vanilla ${1}/03_Combine_and_smooth.R $2 $3 $4 $5 $6 $7 $8

