#!/bin/bash
#!
#! SLURM job script for clust1
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#SBATCH --job-name=22_Plot_ECNF
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0:20:00
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
echo "Markdown script: $1"
echo "Input file: $2"
echo "Output folder: $3"

# Add paths directly in Rmd as I'm too stupid to pass arguments to an Rmd script.
sed -i "s|^THISFILE.*|THISFILE=\"${2}\"|" $1
sed -i "s|^THISPATH.*|THISPATH=\"${3}\"|" $1

Rscript -e "rmarkdown::render('$1')"

# Go back to original file
sed -i "s|^THISFILE.*|THISFILE=\"\"|" $1
sed -i "s|^THISPATH.*|THISPATH=\"\"|" $1
