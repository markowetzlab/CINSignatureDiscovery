#!/bin/bash
#!
#! SLURM job script for clust1
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#SBATCH --job-name=23_Fit_MixtureModels
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
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
echo "Path to main_function script: $2"
echo "Input ECNF file: $3"
echo "Output path: $4"
echo "Feature distribution to fit: $5"
echo "Family of distribution use to fit: $6"
echo "Minimum number of components: $7"
echo "Maximum number of components: $8"
echo "Prior: $9"
echo "Number of iterations: ${10}"
echo "Model decision criterion: ${11}"
echo "Seed for this run: ${12}"

Rscript --vanilla ${1}/23_Fit_MixtureModels.R $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12}



  
