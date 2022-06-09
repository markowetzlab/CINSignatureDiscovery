#!/bin/bash
#!
#! SLURM job script for clust1
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#SBATCH --job-name=31_Calculate_SxC_Matrix
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
echo "Input ECNF: $2"
echo "Input Models: $3"
echo "Output matrix: $4"
echo "Using uninformative priors?: $5"

Rscript --vanilla ${1}/31_Calculate_SxC_Matrix.R $2 $3 $4 $5
