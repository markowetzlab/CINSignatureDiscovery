#!/bin/bash
#!
#! SLURM job script for clust1
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#SBATCH --job-name=11_Compute_OV_sigs_on_TCGA
#SBATCH --time=0:30:00
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
echo "Prepath: $2"
echo "Output folder: $3"
echo "OV signatures: $4"
echo "Geoff' parameters: $5"
echo "Cores: $6"
echo "Geoff's methods: $7"
echo "New methods: $8"
echo "Kerstin's metadata file: $9"
echo "Geoff's extracted features: ${10}"
echo "New filtered data: ${11}"
echo "New smoothed and filtered data: ${12}"
echo "Date: ${13}"

Rscript --vanilla ${1}/11_Compute_OV_sigs_on_TCGA.R $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13}
