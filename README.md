# GSVA Python CLI

Python hooks for R's GSVA bioconductor package to make working in Pandas easier, and a handy CLI for execution of GSVA.

Autodoc manual is here:  https://jason-weirather.github.io/GSVA/

##### Disclaimer

I am not the creator or author of GSVA.  This is a CLI and python hook created to make their package easy to use from the command line and python.

##### This is not the offical site for the GSVA bioconductor package

Find the official R package here

https://doi.org/doi:10.18129/B9.bioc.GSVA

##### And if you find this useful, please cite the author's publication

Hänzelmann S, Castelo R and Guinney J (2013). “GSVA: gene set variation analysis for microarray and RNA-Seq data.” BMC Bioinformatics, 14, pp. 7. doi: 10.1186/1471-2105-14-7, http://www.biomedcentral.com/1471-2105/14/7.

## Get GSVA Python CLI

#### Method 1: Install on your system

1. Install R https://www.r-project.org/ 
2. Install the R bioconductor packaqge GSEABase and GSVA 

```
$ Rscript -e 'source("http://bioconductor.org/biocLite.R");\
              library(BiocInstaller);\
              biocLite(pkgs=c("GSEABase","GSVA"),dep=TRUE)'
```

3. Install this package `$ pip install GSVA`

#### Method 2: Run GSVA via the docker

`$ docker pull vacation/gsva:latest`

## Use GSVA Python CLI in your python code

First install GSVA Python CLI on your system as described above. For details on the `gsva(expression_df,genesets_df,...)` function parameters see https://jason-weirather.github.io/GSVA/ 

```
import pandas as pd
from GSVA import gsva

pathways = gsva(expression_df,genesets_df)

```