#!/bin/bash
#!
#! SLURM job script for clust1
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#SBATCH --job-name=23_Fit_MixtureModels
#SBATCH --time=144:00:00
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
echo "Input ECNF file: $2"
echo "Output path: $3"
echo "Feature distribution to fit: $4"
echo "Minimum number of components: $5"
echo "Maximum number of components: $6"
echo "Number of iterations: $7"
echo "Number of cores: $8"

Rscript --vanilla ${1}/23_Fit_MixtureModels.R $2 $3 $4 $5 $6 $7 $8



  
