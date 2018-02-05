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

### Example convert a per-sample per-gene expression matrix to a per-sample per-pathway enrichment matrix

Consider this Jupyter notebook workflow

```
import pandas as pd
from GSVA import gsva, gmt_to_dataframe
```

Read in a Broad reference pathway gmt file.  Notice the "member" and "name" fields.  If you make your own dataframe to use, these are the required column names.

```
genesets_df = gmt_to_dataframe('c2.cp.v6.0.symbols.gmt')
genesets_df.head()
```

|	| description	                                    | member | name                            |
|---|---------------------------------------------------|--------|---------------------------------|
| 0	| http://www.broadinstitute.org/gsea/msigdb/card... | ACSS2  | KEGG_GLYCOLYSIS_GLUCONEOGENESIS |
| 1	| http://www.broadinstitute.org/gsea/msigdb/card... | GCK    | KEGG_GLYCOLYSIS_GLUCONEOGENESIS |
| 2	| http://www.broadinstitute.org/gsea/msigdb/card... | PGK2   | KEGG_GLYCOLYSIS_GLUCONEOGENESIS |
| 3	| http://www.broadinstitute.org/gsea/msigdb/card... | PGK1   | KEGG_GLYCOLYSIS_GLUCONEOGENESIS |
| 4	| http://www.broadinstitute.org/gsea/msigdb/card... | PDHB   | KEGG_GLYCOLYSIS_GLUCONEOGENESIS |

This example has 200 samples

```
expression_df = pd.read_csv('example_expression.csv',index_col=0)
expression_df.iloc[0:5,0:5]
```

| gene_name | S-1    | S-2    | S-3    | S-4    | S-5    |
|-----------|--------|--------|--------|--------|--------|
| MT-CO1    | 13.852 | 12.328 | 13.055 | 11.898 | 10.234 |
| MT-CO2    | 13.406 | 12.383 | 13.281 | 11.578 | 11.156 |
| MT-CO3    | 13.234 | 12.109 | 13.352 | 11.531 | 10.422 |
| MT-ATP8   | 13.805 | 11.789 | 13.414 | 11.883 | 11.141 |
| MT-ATP6   | 13.500 | 11.703 | 13.227 | 11.219 | 10.836 |

The default command runs without verbose message output. but take notice, that genes that are not part of the `expression_df` are dropped from the analysis, and depending on your choice of GSVA method, genes for which there is not enough expression (i.e. all zero expression) will be dropped.

```
pathways_df = gsva(expression_df,genesets_df)
pathways_df.iloc[0:5,0:5]
```

| name                    | S-1       | S-2       | S-3       | S-4      | S-5       |
|-------------------------|-----------|-----------|-----------|----------|-----------|
| BIOCARTA_41BB_PATHWAY   | 0.068631  | 0.257169  | -0.146907 | 0.020151 | -0.234537 |
| BIOCARTA_ACE2_PATHWAY   | 0.110822  | -0.222310 | -0.161572 | 0.370659 | -0.003318 |
| BIOCARTA_ACH_PATHWAY    | 0.514193  | 0.149291  | 0.226279  | 0.289960 | 0.016071  |
| BIOCARTA_ACTINY_PATHWAY | -0.014494 | 0.407871  | -0.062163 | 0.055607 | 0.424726  |
| BIOCARTA_AGPCR_PATHWAY  | 0.622482  | -0.012845 | 0.317349  | 0.286368 | 0.022540  |
