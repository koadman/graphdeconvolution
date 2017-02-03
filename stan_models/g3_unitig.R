#!/usr/bin/env Rscript
library(rstan)
rstan_options(auto_write=TRUE)

args <- commandArgs(trailingOnly = F)  
script.dir <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
stancode<-paste(script.dir,"genotypes3_unitig_graph.stan",sep="/")

sm <- stan_model(stancode,model_name="g3_unitig")

source(commandArgs(TRUE)[1])
dat <- list(G=G, maxdepth=maxdepth,V=V,S=S,depths=depths,adj1count=adj1count,adj1source=adj1source,adj1dest=adj1dest,adj2count=adj2count,adj2source=adj2source,adj2dest1=adj2dest1,adj2dest2=adj2dest2,lengths=lengths,expected_length=expected_length)

vbfit <- vb(sm, data = dat, tol_rel_obj=0.005, iter=5000, sample_file = 'geno.csv', algorithm="meanfield")

