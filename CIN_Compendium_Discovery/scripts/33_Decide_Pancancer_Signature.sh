#!/bin/bash
#!
#! SLURM job script for clust1
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#SBATCH --job-name=33_Decide_Pancancer_Signature
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
echo "Path to BayesNMF files: $2"
echo "Signature Matrix (W) file pattern: $3"
echo "Exposure Matrix (H) file pattern: $4"
echo "Log file pattern: $5"
echo "Correlation threshold: $6"
echo "Output directory: $7"
echo "SxC as input -> transposed. CxS as input -> untransposed. Transposed: $8"
echo "On what decision criteria should the optimal solution be based on: $9"
echo "Which K to analyse: ${10}"
echo "Path to metadata file: ${11}"
echo "Path to cancer colours: ${12}"

Rscript --vanilla ${1}/33_Decide_Pancancer_Signature.R $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12}
