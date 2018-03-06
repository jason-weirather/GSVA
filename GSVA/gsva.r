#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)
#if (length(args)<15) {
#  stop("Must supply inputs\n",call.=FALSE)
#}

library(GSEABase)
library(GSVA)


method=args[1]
kcdf=args[2]
abs.ranking=as.logical(toupper(args[3]))
min.sz=as.numeric(args[4])
max.sz=args[5]
if(max.sz=='None') {
    max.sz = Inf
} else {
    max.sz = as.numeric(max.sz)
}
parallel.sz=as.numeric(args[6])
parallel.type=args[7]
mx.diff=as.logical(toupper(args[8]))
tau=args[9]
ssgsea.norm=as.logical(toupper(args[10]))
verbose=as.logical(toupper(args[11]))
tempdir=args[12]


geneSets = getGmt(file.path(tempdir,'gs.gmt'))
mat = as.matrix(read.csv(file.path(tempdir,'expr.csv'),header=TRUE,row.names=1,check.names=FALSE))

output = NA
if(tau != 'None') {
    if(tau=='NA') { tau=NA}
    else { tau = as.numeric(tau)}
    output = gsva(mat,geneSets,
        method=method,
        kcdf=kcdf,
        abs.ranking=abs.ranking,
        min.sz=min.sz,
        max.sz=max.sz,
        parallel.sz=parallel.sz,
        parallel.type=parallel.type,
        mx.diff=mx.diff,
        tau=tau,
        ssgsea.norm=ssgsea.norm,
        verbose=verbose
    )
} else {
    output = gsva(mat,geneSets,
        method=method,
        kcdf=kcdf,
        abs.ranking=abs.ranking,
        min.sz=min.sz,
        max.sz=max.sz,
        parallel.sz=parallel.sz,
        parallel.type=parallel.type,
        mx.diff=mx.diff,
        ssgsea.norm=ssgsea.norm,
        verbose=verbose
    )
}
ofile = file.path(tempdir,"pathways.csv")
print(ofile)
write.csv(output,ofile)
