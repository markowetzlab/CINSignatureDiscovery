# Prepare code and environment

## Download SignatureAnalyzer-GPU

BASE="/home/drews01/SignatureDiscovery/Denovo_Signature_Discovery"
cd $BASE

git clone https://github.com/broadinstitute/SignatureAnalyzer-GPU
# We want to checkout the version the original paper used
git checkout 8ce68d08

## Prepare conda environment
conda create --name bayesnmf --file conda_env_specs.txt


## Test
cd SignatureAnalyzer-GPU
conda activate bayesnmf

python SignatureAnalyzer-GPU.py --data example_data/POLEMSI_counts_matrix.txt --max_iter=100000 --output_dir POLEMSI_EXAMPLE --prior_on_W L1 --prior_on_H L2 --labeled

# This should run now. If not, please revisit the conda environment and software versions.
