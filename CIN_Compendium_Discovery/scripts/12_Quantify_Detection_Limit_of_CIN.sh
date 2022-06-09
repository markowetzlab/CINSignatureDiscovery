#!/bin/bash
#!
#! SLURM job script for clust1
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#SBATCH --job-name=12_Quantify_Detection_Limit_of_CIN
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=0:15:00
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
echo "Absolute exposures: $2"
echo "Relative exposures: $3"
echo "OV signatures: $4"
echo "Signature detection limit: $5"
echo "Kerstin's metadata file: $6"
echo "Clinical data: $7"
echo "Output folder: $8"

Rscript --vanilla ${1}/12_Quantify_Detection_Limit_of_CIN.R $2 $3 $4 $5 $6 $7 $8
