# Multi-Omics MIX 
## Benchmark of multi-omics joint Dimensionality Reduction (jDR) approaches in cancer study

We here extensively benchmark 9 representative jDR approaches in three contexts: 
1. samples clustering from multi-omics simulated data
2. ability to identify factors associated with survival or clinical annotations and metagenes associated with biological annotations (Reactome, GO, Hallmarks) in bulk multi-omics TCGA data from 10 cancer types.
3. cells clustering based on scRNA-seq and scATAC-seq data from three cell lines.

The benchmarked methods are:
* iCluster 
* Integrative NMF (intNMF) 
* Joint and individual variation explained (JIVE) 
* Multiple co-inertia analysis (MCIA) 
* Multi-Omics Factor Analysis (MOFA)
* Multi-Study Factor Analysis (MSFA) 
* Regularized Generalized Canonical Correlation Analysis (RGCCA) 
* matrix-tri-factorization (scikit-fusion) 
* tensorial Independent Component Analysis (tICA)

Please note that due to long running time, MSFA is not executed by default.

Each of the three sub-benchmarks above corresponds to a differnet Jupyter notebook in this repositiory:
1. 'Comparison in simulated data.ipynb'
2. 'Comparison in cancer data.ipynb'
3. 'Comparison in single-cell data.ipynb'

## Input data

The input data for the three notebooks are organized as follows: 
1. Comparison in simulated data: data simulated as part of the notebook 
2. Comparison in cancer data: bulk cancer data from TCGA, which was used for a previous study, available at
  http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html. A download script
  is provided to retrieve this dataset
3. Comparison in single-cell data: single-cell data from Liu et al. (Nat Commun. 2019, 10(1):470) are available in `./data/` folder

For running 2. also annotations from MSigDB are required and a download script is provided to retrieve them.

## Install the software environment

* Install conda from https://docs.conda.io/en/latest/miniconda.html
 * create a new environment: `conda create -n momix -c conda-forge -c bioconda -c lcantini momix r-irkernel`
* download annotations and data with the two provided shell scripts (should work on linux and OSX):
 `./download-annotations.sh` and `./download-data.sh`.


## Run the notebooks

* Enter the conda environment: `conda activate momix`.
* Launch the notebook with `jupyter-notebook`.
