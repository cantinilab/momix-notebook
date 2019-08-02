# A mix of multi-omics data

These notebooks are applied on 3 different datasets:

* simulated data, which is generated on demand
* bulk cancer data from TCGA, which was used for a previous study, available at
  http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html. A download script
  is provided to retrieve this dataset
* single-cell data from Liu et al. (Nat Commun. 2019, 10(1):470). The authors
  provide the pre-processed data on-demand.
* annotations from MSigDB. A download script is provided to retrieve them.


## Install the software environment

* Install conda from https://docs.conda.io/en/latest/miniconda.html
 * create a new environment: `conda create -n momix -c conda-forge -c bioconda -c lcantini momix r-irkernel`
* download annotations and data with the two provided shell scripts (should work on linux and OSX):
 `./download-annotations.sh` and `./download-data.sh`.


## Run the notebooks

* Enter the conda environment: `conda activate momix`.
* Launch the notebook with `jupyter-notebook`.
