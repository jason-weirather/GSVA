"""Execute bioconductors GSVA transformation of gene expression into pathway enrichment

"""
import argparse, sys, os
import pandas as pd 
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE

def gsva(df,gmt_df=None,gmt_file=None,
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
    if not tempdir:
        tempdir =  smkdtemp(prefix="weirathe.",dir=gettempdir().rstrip('/'))
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

def cli():
    args = do_inputs()
    # Now read in the input files for purposes of standardizing inputs
    df = None
    if args.tsv_in:
        df = pd.read_csv(args.input,sep="\t",index_col=0)
    else:
        df = pd.read_csv(args.input,index_col=0)
    gmt = gmt_to_pd(args.gmt)
    result = gsva(df,gmt_df=gmt,
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
    if args.output:
        result.to_csv(args.output,args)
    else:
        result.to_csv(os.path.join(args.tempdir,'final.csv'))
        with open(os.path.join(args.tempdir,'final.csv')) as inf:
            for line in inf:
                sys.stdout.write(line)

def gmt_to_pd(fname):
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

def main(args):
  inf = sys.stdin
  of = sys.stdout
  if args.input != '-':
    if args.input[-3:] == '.gz': inf = gzip.open(args.input)
    else: inf = open(args.input)
  if args.output:
    of = open(args.output,'w')
  stream = FASTAStream(inf)
  for fa in stream:
    if re.match('[\t]',fa.header):
      sys.stderr.write("ERROR: tab in header cannot convert to tsv")
      sys.stderr.write("\n")
    of.write(fa.header+"\t"+fa.sequence.replace("\n",'')+"\n")
  of.close()

def do_inputs():
    # Setup command line inputs
    parser=argparse.ArgumentParser(description="Execute R bioconductors GSVA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group0 = parser.add_argument_group('Input options')
    group0.add_argument('input',help="Use - for STDIN")
    group0.add_argument('--tsv_in',action='store_true',help="Exepct CSV by default, this overrides to tab")
    group0.add_argument('--gmt',required=True,help='GMT file with pathways')

    group2 = parser.add_argument_group('Output options')
    group2.add_argument('--csv_out',action='store_true',help="Exepct CSV output format with comma delimmiter and double quoted text")
    group2.add_argument('--output','-o',help="Specifiy path to write transformed data")
    group2.add_argument('--meta_output',help="Speciify path to output additional run information")

    group1 = parser.add_argument_group('command options')
    method_str = '''
Method to employ in the estimation of gene-set enrichment scores per sample.
By default this is set to gsva (Hänzelmann et al, 2013) and other options
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
as ’highly’ activated.
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
both the gsva (Hänzelmann et al., 2013) and the ssgsea (Barbie et al., 2009)
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

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)

