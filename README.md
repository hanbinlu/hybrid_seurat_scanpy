Use Scanpy (Python) and Seurat (R) to analyze single cell RNA-seq data. Seurat is used for data normalization (`SCTransform`), and integration of different stages (`FindIntegrationAnchors` and `IntegrateData`). Scanpy is used for advanced embedding (besides PCA and UMAP), data visualization. 
AnnData (PY) is used as the datastructure to handle the single cell data. To communicate between seurat object and anndata object (of Scanpy), Anndata2ri is used to convert. This is done in python enabled by rpy2 to embed R in python.

To set up the above workflow, we installed Seura, Anndata2ri, Scanpy using Conda.

```
# R=4.0.5 is required for compatibilities between Anndata2ri=1.0.5 and rpy2=3.4.2
# see https://github.com/theislab/anndata2ri/issues/63
conda create -n scRNA -c conda-forge r-essentials r-base=4.0
conda activate scRNA

# install seurat
# start R
install.packages("Seurat")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("glmGamPoi")

# Anndata, Anndata2ri, and Scanpy
pip install 'rpy2==3.4.2'
pip install anndata2ri
conda install seaborn scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph leidenalg 
conda install -c bioconda scanpy

# jupyter
conda install ipykernel nb_black 
python -m ipykernel install --user --name scrna --display-name "scrna"
```