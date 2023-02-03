# Discovery of copy number signatures

This repo will guide you through the process of how Drews et al. (2022) have obtained signatures from copy number features. There are two main folders, one is containing the code and the data detailing the step-by-step process of how the published set of 17 signatures have been discovered (see [**CIN Compendium Discovery**](CIN_Compendium_Discovery/)). The other folder contains the code and a brief summary of how to obtain de-novo signatures from copy number profiles (see [**Denovo Signature Discovery**](Denovo_Signature_Discovery/)).


## CIN Compendium Discovery

The workflow is currently optimised for a slurm HPC environment. 

The main file is `scripts/00_WORKFLOW_SLURM.sh`. This file lists for each step which files are needed and which output are to be expected.


## De-novo Signature Discovery

To obtain signatures from copy number profiles, you will need to extract copy number features, calculate the posterior probability to the 43 feature components and sum the values for each patient / sample. This conversion of the copy number profiles to a 43D vector is a bespoke process but essentially a pre-processing step before the signature discovery can take place. Therefore, its code plus much more documentation can be found in the R package and the github repository [**SignatureQuantification**](https://github.com/markowetzlab/CINSignatureQuantification).

Once you have obtained a matrix, you can use the code in this repository to identify signatures. For this code to run, you will need [SignatureAnalyzer-GPU](https://github.com/broadinstitute/SignatureAnalyzer-GPU) and have Pytorch (v1.4.0, Python 3.6) and CUDA toolkit installed (v10.0.130, Python 3.6). The easiest way is to create a conda environment called **bayesnmf** based on this specifications file: 
```
conda create --name bayesnmf --file Denovo_Signature_Discovery/conda_env_specs.txt
```

### Example data

You will need to have `git` and `conda` installed for the code to be working. If so, run this code to create the conda environment and get the code for SignatureAnalyzer-GPU:

```
./Denovo_Signature_Discovery/scripts/0_Prepare_code_and_env.sh
```

Then call SignatureAnalyzer-GPU on the example data:
```
./Denovo_Signature_Discovery/scripts/1_Signature_derivation_from_input_matrix.sh`
```

If you want to use it for your own data, open the file `Denovo_Signature_Discovery/scripts/1_Signature_derivation_from_input_matrix.sh` and edit the paths pointing towards the input file and output directory.


Then call the following R script from the command line to identify and extract the optimal solution:
```
Rscript Denovo_Signature_Discovery/scripts/2_Decide_signature_solution.R
```

When running this code on your own data you will want to change the variables `PATHTOFILES`, `OUTPUTDIR` and `OVERRIDEK` to match your settings.

## Maintenance

For any issues please contact please contact Florian Markowetz Florian.Markowetz@cruk.cam.ac.uk and Geoff Macintyre gmacintyre@cnio.es.


## Contact

If you experience any issues or have questions about the code, please open a Github issue with a minimum reproducible example. For questions around collaborations or sensitive patient data, please contact us directly at Florian Markowetz <Florian.Markowetz@cruk.cam.ac.uk> and Geoff Macintyre <gmacintyre@cnio.es>.

## Licence
The contents of this repository are copyright (c) 2022, University of Cambridge and Spanish National Cancer Research Centre (CNIO).

The contents of this repository are published and distributed under the GAP Available Source License v1.0 (ASL). 

The contents of this repository are distributed in the hope that it will be useful for non-commercial academic research, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the ASL for more details. 

The methods implemented in the code are the subject of pending patent application GB 2114203.9.

Any commercial use of this code is prohibited.

