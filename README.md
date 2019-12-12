# Multi-Omics MIX (MOMIX)
## Benchmark of multi-omics joint Dimensionality Reduction (jDR) approaches in cancer study

![image](https://drive.google.com/uc?export=view&id=11MC0sJ3rPvqZcoaBfudZ5X-A32OPU1iM)

We here extensively benchmark 9 representative jDR approaches in three contexts: 
1. samples clustering from multi-omics simulated data
2. ability to identify factors associated with survival or clinical annotations and metagenes associated with biological annotations (Reactome, GO, Hallmarks) in bulk multi-omics TCGA data from 10 cancer types.
3. cells clustering based on scRNA-seq and scATAC-seq data from three cell lines.

The benchmarked methods are:
* [iCluster](https://cran.r-project.org/web/packages/iCluster/index.html)
* [Integrative NMF (intNMF)](https://cran.r-project.org/web/packages/IntNMF/index.html) 
* [Joint and individual variation explained (JIVE)](https://cran.r-project.org/web/packages/r.jive/index.html) 
* [Multiple co-inertia analysis (MCIA)](https://bioconductor.org/packages/release/bioc/html/omicade4.html) 
* [Multi-Omics Factor Analysis (MOFA)](https://github.com/bioFAM/MOFA)
* [Multi-Study Factor Analysis (MSFA)](https://github.com/rdevito/MSFA) 
* [Regularized Generalized Canonical Correlation Analysis (RGCCA)](https://cran.r-project.org/web/packages/RGCCA/index.html) 
* [matrix-tri-factorization (scikit-fusion)](https://github.com/marinkaz/scikit-fusion) 
* [tensorial Independent Component Analysis (tICA)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1455-8)

Please note that due to long running time, MSFA is not executed by default.

Each of the three sub-benchmarks above corresponds to a differnet Jupyter notebook in this repositiory:
1. `Comparison in simulated data.ipynb`
2. `Comparison in cancer data.ipynb`
3. `Comparison in single-cell data.ipynb`

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

## Run the notebooks

* Enter the conda environment: `conda activate momix`.
* Launch the notebook with `jupyter-notebook`.

## Compare the jDR methods on new data
* If the inputs are bulk multi-omics data then run the 'Comparison in cancer data.ipynb' notebook changing the folder where to access the input data from `data/cancer` to `data/folder_new_data`.
* If the inputs are single-cell multi-omics data  then run the 'Comparison in single-cell data.ipynb' notebook changing the folder where to access the input data from `data/single-cell` to `data/folder_new_data`. 

##  Compare a new jDR method in respect to the state-of-the-art
Add the commans to run the new algorithm in the end of the function `scripts/ranfactorization.R`. In this case the output of that function will have to contain the name of the new added method, the factors obtained with the method and the list of weight matrices obtained with the new method.
