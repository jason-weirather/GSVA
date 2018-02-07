# GSVA / ssGSEA command-line interface and Python module

The GSVA (gene-set variance analysis) package from R bioconductor provides efficient computation of single-sample gene-set enrichment analysis (ssGSEA). This pakcage provides a python implmented CLI, and Python module with Pandas inputs and outputs, as well as a docker to run this R package.

* Repository is here: https://github.com/jason-weirather/GSVA
* Autodoc manual is here:  https://jason-weirather.github.io/GSVA/

##### Disclaimer

I am not the creator or author of GSVA.  This is a CLI and python hook created to make their package easy to use from the command line and python. *This is not the offical site for the GSVA bioconductor package*

Find the official R package here

https://doi.org/doi:10.18129/B9.bioc.GSVA

##### And if you find this useful, please cite the author's publication

Hänzelmann S, Castelo R and Guinney J (2013). “GSVA: gene set variation analysis for microarray and RNA-Seq data.” BMC Bioinformatics, 14, pp. 7. doi: 10.1186/1471-2105-14-7, http://www.biomedcentral.com/1471-2105/14/7.

## Quickstart - CLI through docker

### Execute GSVA in docker

Just be careful to let docker know all the volumes you need to mount.  Here we will do everything in our current directory.

1. Start with an expression csv with gene-wise rows and sample-wise columns

```
$ cat example_expression.csv | cut -f 1-3 -d ',' | head -n 6 
gene_name,S-1,S-2
MT-CO1,13.852,12.328
MT-CO2,13.405999999999999,12.383
MT-CO3,13.234000000000002,12.109000000000002
MT-ATP8,13.805,11.789000000000001
MT-ATP6,13.5,11.703
```

2. Use any gene sets in **gmt** format where each row follows the convention `name <tab> description <tab> gene1 <tab> gene2 <tab> ... <tab> geneN`

```
cat c2.cp.v6.0.symbols.gmt | head -n 6 | cut -f 1-5
KEGG_GLYCOLYSIS_GLUCONEOGENESIS	http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_GLYCOLYSIS_GLUCONEOGENESIS	ACSS2	GCK	PGK2
KEGG_CITRATE_CYCLE_TCA_CYCLE	http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_CITRATE_CYCLE_TCA_CYCLE	IDH3B	DLST	PCK2
KEGG_PENTOSE_PHOSPHATE_PATHWAY	http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_PENTOSE_PHOSPHATE_PATHWAY	RPE	RPIA	PGM2
KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS	http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS	UGT1A10	UGT1A8	RPE
KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM	http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM	MPI	PMM2	PMM1
KEGG_GALACTOSE_METABOLISM	http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_GALACTOSE_METABOLISM	GCK	GALK1	GLB1
```

3. Run GSVA

```
$ docker run -v $(pwd):$(pwd) vacation/gsva:1.0.4 \
    GSVA --gmt $(pwd)/c2.cp.v6.0.symbols.gmt \
         $(pwd)/example_expression.csv \
         --output $(pwd)/example_pathways.csv
```

##### You're done.  Thats it.  Enjoy, check your output.

Running outside of docker on your system is just as easy (actually easier) but you need to have the required programs installed (see below). 

```
$ cat example_pathways.csv | cut -f 1-3 -d ',' | head -n 6
name,S-1,S-2
BIOCARTA_41BB_PATHWAY,0.0686308398590944,0.257169127694153
BIOCARTA_ACE2_PATHWAY,0.11082238459933501,-0.22231034473486602
BIOCARTA_ACH_PATHWAY,0.514192767265737,0.149291024991685
BIOCARTA_ACTINY_PATHWAY,-0.0144944990305252,0.407870971955071
BIOCARTA_AGPCR_PATHWAY,0.6224821629523449,-0.0128449355033173
```


### For more advanced options you can list the options

```
$ docker run vacation/gsva:1.0.4 GSVA -h
usage: GSVA [-h] [--tsv_in] --gmt GMT [--tsv_out] [--output OUTPUT]
            [--meta_output META_OUTPUT] [--method {gsva,ssgsea,zscore,plage}]
            [--kcdf {Gaussian,Poisson,none}] [--abs_ranking] [--min_sz MIN_SZ]
            [--max_sz MAX_SZ] [--parallel_sz PARALLEL_SZ]
            [--parallel_type PARALLEL_TYPE] [--mx_diff MX_DIFF] [--tau TAU]
            [--ssgsea_norm SSGSEA_NORM] [--verbose]
            [--tempdir TEMPDIR | --specific_tempdir SPECIFIC_TEMPDIR]
            input

Execute R bioconductors GSVA

optional arguments:
  -h, --help            show this help message and exit

Input options:
  input                 Use - for STDIN
  --tsv_in              Exepct CSV by default, this overrides to tab (default:
                        False)
  --gmt GMT             GMT file with pathways (default: None)

Output options:
  --tsv_out             Override the default CSV and output TSV (default:
                        False)
  --output OUTPUT, -o OUTPUT
                        Specifiy path to write transformed data (default:
                        None)
  --meta_output META_OUTPUT
                        Speciify path to output additional run information
                        (default: None)

command options:
  --method {gsva,ssgsea,zscore,plage}
                        Method to employ in the estimation of gene-set
                        enrichment scores per sample. By default this is set
                        to gsva (Hanzelmann et al, 2013) and other options 6
                        gsva are ssgsea (Barbie et al, 2009), zscore (Lee et
                        al, 2008) or plage (Tomfohr et al, 2005). The latter
                        two standardize first expression profiles into
                        z-scores over the samples and, in the case of zscore,
                        it combines them together as their sum divided by the
                        square-root of the size of the gene set, while in the
                        case of plage they are used to calculate the singular
                        value decomposition (SVD) over the genes in the gene
                        set and use the coefficients of the first right-
                        singular vector as pathway activity profile. (default:
                        gsva)
  --kcdf {Gaussian,Poisson,none}
                        Character string denoting the kernel to use during the
                        non-parametric estimation of the cumulative
                        distribution function of expression levels across
                        samples when method="gsva". By default,
                        kcdf="Gaussian" which is suitable when input
                        expression values are continuous, such as microarray
                        fluorescent units in logarithmic scale, RNA-seq log-
                        CPMs, log-RPKMs or log-TPMs. When input expression
                        values are integer counts, such as those derived from
                        RNA-seq experiments, then this argument should be set
                        to kcdf="Poisson". This argument supersedes arguments
                        rnaseq and kernel, which are deprecated and will be
                        removed in the next release. (default: Gaussian)
  --abs_ranking         Flag used only when mx_diff=TRUE. When
                        abs_ranking=FALSE (default) a modified Kuiper
                        statistic is used to calculate enrichment scores,
                        taking the magnitude difference between the largest
                        positive and negative random walk deviations. When
                        abs.ranking=TRUE the original Kuiper statistic that
                        sums the largest positive and negative random walk
                        deviations, is used. In this latter case, gene sets
                        with genes enriched on either extreme (high or low)
                        will be regarded as'highly' activated. (default:
                        False)
  --min_sz MIN_SZ       Minimum size of the resulting gene sets. (default: 1)
  --max_sz MAX_SZ       Maximum size of the resulting gene sets. (default:
                        None)
  --parallel_sz PARALLEL_SZ
                        Number of processors to use when doing the
                        calculations in parallel. This requires to previously
                        load either the parallel or the snow library. If
                        parallel is loaded and this argument is left with its
                        default value (parallel_sz=0) then it will use all
                        available core processors unless we set this argument
                        with a smaller number. If snow is loaded then we must
                        set this argument to a positive integer number that
                        specifies the number of processors to employ in the
                        parallel calculation. (default: 0)
  --parallel_type PARALLEL_TYPE
                        Type of cluster architecture when using snow.
                        (default: SOCK)
  --mx_diff MX_DIFF     Offers two approaches to calculate the enrichment
                        statistic (ES) from the KS random walk statistic.
                        mx_diff=FALSE: ES is calculated as the maximum
                        distance of the random walk from 0. mx_diff=TRUE
                        (default): ES is calculated as the magnitude
                        difference between the largest positive and negative
                        random walk deviations. (default: True)
  --tau TAU             Exponent defining the weight of the tail in the random
                        walk performed by both the gsva (Hanzelmann et al.,
                        2013) and the ssgsea (Barbie et al., 2009) methods. By
                        default, this tau=1 when method="gsva" and tau=0.25
                        when method="ssgsea" just as specified by Barbie et
                        al. (2009) where this parameter is called alpha.
                        (default: None)
  --ssgsea_norm SSGSEA_NORM
                        Logical, set to TRUE (default) with method="ssgsea"
                        runs the SSGSEA method from Barbie et al. (2009)
                        normalizing the scores by the absolute difference
                        between the minimum and the maximum, as described in
                        their paper. When ssgsea_norm=FALSE this last
                        normalization step is skipped. (default: True)
  --verbose             Gives information about each calculation step.
                        (default: False)

Temporary folder parameters:
  --tempdir TEMPDIR     The temporary directory is made and destroyed here.
                        (default: /tmp)
  --specific_tempdir SPECIFIC_TEMPDIR
                        This temporary directory will be used, but will remain
                        after executing. (default: None)
```

## Installation

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

### Workflow example - Go from an expression-based tSNE plot to a pathway-based tSNE plot in a Jupyter notebook

Here we will convert a per-sample per-gene expression matrix to a per-sample per-pathway enrichment matrix. We will plot the values using tSNE.

These code snipits and outputs are from a Jupyter notebook.


```python
import pandas as pd
from GSVA import gsva, gmt_to_dataframe
# Some extras to look at the high dimensional data
from plotnine import *
from sklearn.manifold import TSNE
```

Read in a Broad reference pathway gmt file.  Notice the "member" and "name" fields.  If you make your own dataframe to use, these are the required column names.

```python
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

```python
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

```python
XV = TSNE(n_components=2).\
    fit_transform(expression_df.T)
df = pd.DataFrame(XV).rename(columns={0:'x',1:'y'})
(ggplot(df,aes(x='x',y='y'))
 + geom_point(alpha=0.2)
)
```

![Gene Expression](https://i.imgur.com/Qbwds5H.png)

The default command runs without verbose message output. but take notice, that genes that are not part of the `expression_df` are dropped from the analysis, and depending on your choice of GSVA method, genes for which there is not enough expression (i.e. all zero expression) will be dropped.

```python
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

```python
YV = TSNE(n_components=2).\
    fit_transform(pathways_df.T)
pf = pd.DataFrame(YV).rename(columns={0:'x',1:'y'})
(ggplot(pf,aes(x='x',y='y'))
 + geom_point(alpha=0.2)
)
```

![Pathway Enrichment](https://i.imgur.com/2pxjoRr.png)
