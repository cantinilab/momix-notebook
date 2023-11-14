#initialize conda environment
conda create -n momix -c conda-forge -y  "python>=3.6,<3.7" r-base

# activate env
conda activate momix

# install jupyter through pip
pip install jupyter

# install r dependencies (this will take a while) 
conda install -c conda-forge -y r-survival r-rgcca r-icluster r-gparotation r-ggplot2 r-clustercrit r-devtools r-irkernel

# install missing python dependencies directly from github 
pip install git+https://github.com/mims-harvard/scikit-fusion

# install omicade4
conda install -c bioconda -y bioconductor-omicade4

# install mofa
pip install mofapy
Rscript -e 'devtools::install_github("bioFAM/MOFA")'

#install fgsea from devtoools
Rscript -e 'devtools::install_github("ctlab/fgsea")'

# install MSFA
Rscript -e 'devtools::install_github("lcan88/MSFA")'

#install missing packages frm CRAN:
Rscript -e 'install.packages(c("r.jive","IntNMF", "InterSIM", "tensorBSS"),  repos="http://cran.us.r-project.org")'


#download notebooks
git clone https://github.com/cantinilab/momix-notebook

# run jupyter
jupyter notebook

