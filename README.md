# scAnalyzer 

This package provides an integrative analysis pipeline for single-cell RNA-seq data. It implements methods to perform quality control (QC), doublets removal, normalization, clustering, cell type annotation, malignant cell identification, immature cancer cell identification and data visualization.


## Installation
Installing the package in a fresh R environment may take a long time.

### Prerequisites
To install scAnalyzer, you have to first install conda following the document (https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html#install-macos-silent). Afterwards you have to register the bioconda and conda-forge channels as a package source for conda.

~~~
$ conda config --add channels conda-forge
$ conda config --add channels bioconda
$ conda config --add channels r
~~~

Install R, and the JACS package which is required for infercnv.

~~~
$ conda install -c conda-forge r-base=4.1
$ conda install -c conda-forge r-rjags
~~~

### Installation using R
~~~
> install.packages("BiocManager", repos = "https://cloud.r-project.org")
> BiocManager::install(c("Seurat", "ggplot2", "dplyr", "infercnv", "SingleR", "celldex", "scmap", "navinlabcode/copykat", "chris-mcginnis-ucsf/DoubletFinder", "bm2-lab/scLearn"))
> BiocManager::install("WubingZhang/scAnalyzer")
~~~

## Documentation
~~~
> ?scAnalyzer
> ?CellTypeAnnotator
> ?cnvInfer
~~~


