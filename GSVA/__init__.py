"""Execute bioconductors GSVA transformation of gene expression into pathway enrichment.

This python package gives both a CLI interface and a python module to work with GSVA in Python Pandas DataFrames.

Find the official R package here:

https://doi.org/doi:10.18129/B9.bioc.GSVA

And if you find this useful, cite the authors publication:

Hänzelmann S, Castelo R and Guinney J (2013). "GSVA: gene set variation analysis for microarray and RNA-Seq data.” BMC Bioinformatics, 14, pp. 7. doi: 10.1186/1471-2105-14-7, http://www.biomedcentral.com/1471-2105/14/7.

"""
import argparse, sys, os
import pandas as pd 
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE

def gsva(expression_df,geneset_df=None,
         method='gsva',
         kcdf='Gaussian',
         abs_ranking=False,
         min_sz=1,
         max_sz=None,
         parallel_sz=0,
         parallel_type="SOCK",
         mx_diff=True,
         tau=None,
         ssgsea_norm=True,
         verbose=False,
         tempdir= None
         ):
    """GSVA function for use with pandas DataFrame objects

    :param expression_df: REQUIRED: Expression data indexed on gene names column labels as sample ids
    :type expression_df: pandas.DataFrame
    :param geneset_df: REQUIRED: Genesets and their members in a dataframe
    :type geneset_df: pandas.DataFrame
    :param method: Method to employ in the estimation of gene-set enrichment scores per sample. By default this is set to gsva (Hänzelmann et al, 2013) and other options 6 gsva are ssgsea (Barbie et al, 2009), zscore (Lee et al, 2008) or plage (Tomfohr et al, 2005). The latter two standardize first expression profiles into z-scores over the samples and, in the case of zscore, it combines them together as their sum divided by the square-root of the size of the gene set, while in the case of plage they are used to calculate the singular value decomposition (SVD) over the genes in the gene set and use the coefficients of the first right-singular vector as pathway activity profile.
    :type method: string Default: 'gsva'   
    :param kcdf: Character string denoting the kernel to use during the non-parametric estimation of the cumulative distribution function of expression levels across samples when method="gsva". By default, kcdf="Gaussian" which is suitable when input expression values are continuous, such as microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input expression values are integer counts, such as those derived from RNA-seq experiments, then this argument should be set to kcdf="Poisson". This argument supersedes arguments rnaseq and kernel, which are deprecated and will be removed in the next release.
    :type kcdf: string Default: 'Gaussian'
    :param abs_ranking: Flag used only when mx_diff=TRUE. When abs_ranking=FALSE [default] a modified Kuiper statistic is used to calculate enrichment scores, taking the magnitude difference between the largest positive and negative random walk deviations. When abs.ranking=TRUE the original Kuiper statistic that sums the largest positive and negative random walk deviations, is used. In this latter case, gene sets with genes enriched on either extreme (high or low) will be regarded as 'highly’ activated.
    :type abs_ranking: bool Default: False
    :param min_sz: Minimum size of the resulting gene sets.
    :type min_sz: int Default: 1
    :param max_sz: Maximum size of the resulting gene sets. Leave unset for no limit.
    :type max_sz: int Default: Inf
    :param parallel_sz: Number of processors to use when doing the calculations in parallel. This requires to previously load either the parallel or the snow library. If parallel is loaded and this argument is left with its default value (parallel_sz=0) then it will use all available core processors unless we set this argument with a smaller number. If snow is loaded then we must set this argument to a positive integer number that specifies the number of processors to employ in the parallel calculation.
    :type parallel_sz: int Default: 0
    :param parallel_type: Type of cluster architecture when using snow.
    :type parallel_type: string Default: "SOCK"   
    :param mx_diff: Offers two approaches to calculate the enrichment statistic (ES) from the KS random walk statistic. mx_diff=FALSE: ES is calculated as the maximum distance of the random walk from 0. mx_diff=TRUE (default): ES is calculated as the magnitude difference between the largest positive and negative random walk deviations.
    :type mx_diff: bool Default: True    
    :param tau: Exponent defining the weight of the tail in the random walk performed by both the gsva (Hänzelmann et al., 2013) and the ssgsea (Barbie et al., 2009) methods. By default, this tau=1 when method="gsva" and tau=0.25 when method="ssgsea" just as specified by Barbie et al. (2009) where this parameter is called alpha. Leave unset for defaults.
    :type tau: float    
    :param ssgsea_norm: Logical, set to TRUE (default) with method="ssgsea" runs the SSGSEA method from Barbie et al. (2009) normalizing the scores by the absolute difference between the minimum and the maximum, as described in their paper. When ssgsea_norm=FALSE this last normalization step is skipped.
    :type ssgsea_norm: bool Default: True    
    :param verbose: Gives information about each calculation step.
    :type verbose: bool Default: False
    :param tempdir: Location to write temporary files
    :type tempdir: string Default: System Default
    :returns: pandas.DataFrame
    """
    df = expression_df
    gmt_df = geneset_df

    if not tempdir:
        tempdir =  mkdtemp(prefix="weirathe.",dir=gettempdir().rstrip('/'))
    if verbose:
        sys.stderr.write("Caching to "+tempdir+"\n")
    # Remove genes from the genesets that do not occur in the dataset
    members = gmt_df['member'].unique()
    missing = set(members)-set(df.index)
    original = df.index
    if len(missing) > 0:
        if verbose: sys.stderr.write("WARNING removing "+str(len(missing))+\
          " genes from gene sets that don't exist in the data\n"+\
          ",".join(sorted(list(missing)))+"\n")
    gmt_df = gmt_df[~gmt_df['member'].isin(list(missing))]
    # Write our gene sets
    gmt_df = gmt_df.groupby(['name']).\
        apply(lambda x: "\t".join(sorted(list(x['member'])))).reset_index().rename(columns={0:'members'})
    of = open(os.path.join(tempdir,"gs.gmt"),'w')
    for row in gmt_df.itertuples():
        name = row.name
        description = 'description'
        fields = row.members
        of.write(name+"\t"+description+"\t"+fields+"\n")
    of.close()
    df.to_csv(os.path.join(tempdir,"expr.csv"))
    cur = os.path.dirname(os.path.realpath(__file__))
    rscript = os.path.join(cur,"gsva.r")
    cmd = ["Rscript",rscript]+[str(x) for x in [
    method,kcdf,abs_ranking,min_sz,max_sz,parallel_sz,parallel_type,
    mx_diff,tau,ssgsea_norm,verbose,tempdir
      ]]
    if verbose: sys.stderr.write(" ".join(cmd)+"\n")
    destination = PIPE
    if verbose: destination = sys.stderr
    sp = Popen(cmd,stdout=PIPE,stderr=destination)
    sp.communicate()
    output = pd.read_csv(os.path.join(tempdir,"pathways.csv"),index_col=0)
    output.index.name = 'name'
    return output

def __cli():
    args = __do_inputs()
    # Now read in the input files for purposes of standardizing inputs
    df = None
    if args.tsv_in:
        df = pd.read_csv(args.input,sep="\t",index_col=0)
    else:
        df = pd.read_csv(args.input,index_col=0)
    gmt = gmt_to_dataframe(args.gmt)
    result = gsva(df,geneset_df=gmt,
                  method=args.method,
                  kcdf=args.kcdf,
                  abs_ranking=args.abs_ranking,
                  min_sz=args.min_sz,
                  max_sz=args.max_sz,
                  parallel_sz=args.parallel_sz,
                  parallel_type=args.parallel_type,
                  mx_diff=args.mx_diff,
                  tau=args.tau,
                  ssgsea_norm=args.ssgsea_norm,
                  verbose=args.verbose,
                  tempdir=args.tempdir
                 )
    sep = ','
    if args.tsv_out: sep = "\t"
    if args.output:
        result.to_csv(args.output,sep=sep)
    else:
        result.to_csv(os.path.join(args.tempdir,'final.csv'),sep=sep)
        with open(os.path.join(args.tempdir,'final.csv')) as inf:
            for line in inf:
                sys.stdout.write(line)

def gmt_to_dataframe(fname):
    """ A function to convert gmt files to a pandas dataframe 

    :param fname: path to gmt file
    :type fname: string 
    :returns: pandas.DataFrame

    """
    res = []
    with open(fname) as inf:
        for line in inf:
            f = line.rstrip().split("\t")
            name = f[0]
            description = f[1]
            members = f[2:]
            for member in members:
                res.append(pd.Series({'name':name,'description':description,'member':member}))
    return pd.DataFrame(res)


def __do_inputs():
    # Setup command line inputs
    parser=argparse.ArgumentParser(description="Execute R bioconductors GSVA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group0 = parser.add_argument_group('Input options')
    group0.add_argument('input',help="Use - for STDIN")
    group0.add_argument('--tsv_in',action='store_true',help="Exepct CSV by default, this overrides to tab")
    group0.add_argument('--gmt',required=True,help='GMT file with pathways')

    group2 = parser.add_argument_group('Output options')
    group2.add_argument('--tsv_out',action='store_true',help="Override the default CSV and output TSV")
    group2.add_argument('--output','-o',help="Specifiy path to write transformed data")
    group2.add_argument('--meta_output',help="Speciify path to output additional run information")

    group1 = parser.add_argument_group('command options')
    method_str = '''
Method to employ in the estimation of gene-set enrichment scores per sample.
By default this is set to gsva (Hanzelmann et al, 2013) and other options
6 gsva
are ssgsea (Barbie et al, 2009), zscore (Lee et al, 2008) or plage (Tomfohr
et al, 2005). The latter two standardize first expression profiles into z-scores
over the samples and, in the case of zscore, it combines them together as their
sum divided by the square-root of the size of the gene set, while in the case of
plage they are used to calculate the singular value decomposition (SVD) over
the genes in the gene set and use the coefficients of the first right-singular vector
as pathway activity profile.
    '''
    group1.add_argument('--method',choices=["gsva", "ssgsea", "zscore", "plage"],default='gsva',help=method_str)
    kcdf_str = '''
Character string denoting the kernel to use during the non-parametric estimation
of the cumulative distribution function of expression levels across samples
when method="gsva". By default, kcdf="Gaussian" which is suitable when
input expression values are continuous, such as microarray fluorescent units in
logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input
expression values are integer counts, such as those derived from RNA-seq experiments,
then this argument should be set to kcdf="Poisson". This argument
supersedes arguments rnaseq and kernel, which are deprecated and will be
removed in the next release.
    '''
    group1.add_argument('--kcdf',choices=["Gaussian", "Poisson", "none"],default='Gaussian',help=kcdf_str)
    abs_ranking_str = '''
Flag used only when mx_diff=TRUE. When abs_ranking=FALSE (default) a
modified Kuiper statistic is used to calculate enrichment scores, taking the magnitude
difference between the largest positive and negative random walk deviations.
When abs.ranking=TRUE the original Kuiper statistic that sums the
largest positive and negative random walk deviations, is used. In this latter case,
gene sets with genes enriched on either extreme (high or low) will be regarded
as'highly' activated.
    '''
    group1.add_argument('--abs_ranking',action='store_true',help=abs_ranking_str)
    min_sz_str = '''
Minimum size of the resulting gene sets.
    '''
    group1.add_argument('--min_sz',type=int,default=1,help=min_sz_str)
    max_sz_str = '''
Maximum size of the resulting gene sets.
    '''
    group1.add_argument('--max_sz',type=int,help=max_sz_str)
    parallel_sz_str = '''
Number of processors to use when doing the calculations in parallel. This requires
to previously load either the parallel or the snow library. If parallel is
loaded and this argument is left with its default value (parallel_sz=0) then it
will use all available core processors unless we set this argument with a smaller
number. If snow is loaded then we must set this argument to a positive integer
number that specifies the number of processors to employ in the parallel
calculation.
    '''
    group1.add_argument('--parallel_sz',type=int,default=0,help=parallel_sz_str)
    parallel_type_str = '''
Type of cluster architecture when using snow.
    '''
    group1.add_argument('--parallel_type',default="SOCK",help=parallel_type_str)
    mx_diff_str = '''
Offers two approaches to calculate the enrichment statistic (ES) from the KS
random walk statistic. mx_diff=FALSE: ES is calculated as the maximum distance
of the random walk from 0. mx_diff=TRUE (default): ES is calculated as
the magnitude difference between the largest positive and negative random walk
deviations.
    '''
    group1.add_argument('--mx_diff',type=bool,default=True,help=mx_diff_str)
    tau_str = '''
Exponent defining the weight of the tail in the random walk performed by
both the gsva (Hanzelmann et al., 2013) and the ssgsea (Barbie et al., 2009)
methods. By default, this tau=1 when method="gsva" and tau=0.25 when
method="ssgsea" just as specified by Barbie et al. (2009) where this parameter
is called alpha.
    '''
    group1.add_argument('--tau',type=float,help=tau_str)
    ssgsea_norm_str = '''
Logical, set to TRUE (default) with method="ssgsea" runs the SSGSEA method
from Barbie et al. (2009) normalizing the scores by the absolute difference
between the minimum and the maximum, as described in their paper. When
ssgsea_norm=FALSE this last normalization step is skipped.
    '''
    group1.add_argument('--ssgsea_norm',type=bool,default=True,help=ssgsea_norm_str)
    verbose_str = '''
Gives information about each calculation step.
    '''
    group1.add_argument('--verbose',action='store_true',help=verbose_str)

    # Temporary working directory step 1 of 3 - Definition
    label4 = parser.add_argument_group(title="Temporary folder parameters")
    group3 = label4.add_mutually_exclusive_group()
    group3.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
    group3.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")


    args = parser.parse_args()
    setup_tempdir(args)
    return args  

def setup_tempdir(args):
  if args.specific_tempdir:
    if not os.path.exists(args.specific_tempdir):
      os.makedirs(args.specific_tempdir.rstrip('/'))
    args.tempdir = args.specific_tempdir.rstrip('/')
    if not os.path.exists(args.specific_tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  else:
    args.tempdir = mkdtemp(prefix="weirathe.",dir=args.tempdir.rstrip('/'))
    if not os.path.exists(args.tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  if not os.path.exists(args.tempdir):
    sys.stderr.write("ERROR: Problem creating temporary directory\n")
    sys.exit()
  return

if __name__=="__main__":
  _cli()

