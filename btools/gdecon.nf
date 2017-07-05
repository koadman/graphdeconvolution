#!/usr/bin/env nextflow
/**
 * Usage: gdecon.nf --readlist1=/path/to/pe_readfile_list.txt  --readlist2=/path/to/pe_readfile_list.txt --r1suffix=.r1.fastq.gz
 * (c) 2017 Aaron Darling
 */

readlist1 = Channel.from(file(params.readlist1))
readlist2 = Channel.from(file(params.readlist2))

params.k = 51
params.maxstrains = 25
params.relevance_threshold = 0.05
params.maxtiplen = 100
params.asmcores = 8
params.covk = 31

/*
process clean {
    input:
    file(r) from rfiles1
    output:
    file('*.clean.fq') into readlist

    """
    ${GDECONHOME}/external/bbmap/bbduk.sh in=${r} out=${r}.clean.fq ref=${GDECONHOME}/external/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 qtrim=r trimq=10 tpe tbo
    """
}
*/

// constructs the unitig graph
process unitig {
    input:
    file('rlist1') from readlist1
    file('rlist2') from readlist2
    output:
    file('bothreads.unitigs.fa') into unitigs
    file('bothreads.unitigs.fa') into unitigs2

    """
    cat rlist1 rlist2 >> bothreads.txt
    ${GDECONHOME}/external/gatb/bcalm -in bothreads.txt -kmer-size ${params.k} -abundance-min 2 -nb-cores ${params.asmcores}
    ${GDECONHOME}/external/gatb/btrim bothreads.unitigs.fa ${params.k} ${params.maxtiplen} ${params.asmcores}
    """
}

rfiles1 = Channel.from(file(params.readlist1)).splitText() { it.trim() }
rfiles2 = Channel.from(file(params.readlist2)).splitText() { it.trim() }
rfiles = rfiles1.merge(rfiles2) {o, e -> [o, e]}.spread(unitigs)

// estimates depth of coverage on unitigs
process unitig_cov {
    input:
    set reads,reads2,file('unitigs') from rfiles
    output:
    file('*.counts') into covcounts

    """
    name=`basename $reads ${params.r1suffix}`
    cat $reads $reads2 > bothreads
    ${GDECONHOME}/external/gatb/tigops coverage -out \$name.counts -kmer-size ${params.covk} -reads bothreads -tigs unitigs -name \$name
    """
}

// collect coverage from each time point to a single file
covfefe = covcounts.collect()
process extract_cov {
    input:
    file('*') from covfefe
    output:
    file('coverage.csv') into coverage

"""
#!/usr/bin/env python
import glob
import re
listing = glob.glob('*.counts')
covout = open('coverage.csv', 'w')
tigs = {}
for filename in listing:
    ti = 0
    f = open(filename)
    for line in f:
        if line[0] != '>':
            continue
        m=re.search('mean_(\\d*\\.*\\d*)_', line)
        if not ti in tigs: tigs[ti]=str(ti)
        tigs[ti] += "," + m.group(1)
        ti+=1
for i in range(ti):
    covout.write(tigs[i]+"\\n")

"""
}

// convert the bcalm unitigs to GFA format
process convert_gfa {
    input:
    file('coverage.csv') from coverage
    file('unitigs.fa') from unitigs2

    output:
    file('stan.kmer') into stankmer
    file('stan.kmer') into stankmer2
    file('unitigs.gfa') into gfa

"""
    ${GDECONHOME}/external/gatb/convertToGFA.py unitigs.fa unitigs.gfa ${params.k}
    perl -p -i -e "s/\\s+k:i:/\\tkk:i:/g" unitigs.gfa
    perl -p -i -e "s/KM:f:/km:f:/g" unitigs.gfa
    ${GDECONHOME}/btools/gfa2gdecon.py unitigs.gfa coverage.csv ${params.maxstrains} 1500 > stan.kmer
"""
}

// run the Bayesian NMF in stan
process stankmer_nmf {
    input:
    file('stan.kmer') from stankmer
    each repid from 1,2,3,4
    output:
    file("stan.kmer.*.diag") into stankmerdiag
    file("stan.kmer.*.out") into stankmerout

    """
    ${GDECONHOME}/stan_models/genotypes3_unitig_graph variational output_samples=100 tol_rel_obj=0.008 iter=10000 data file=stan.kmer output file=stan.kmer.${repid}.out diagnostic_file=stan.kmer.${repid}.diag
    """
}


// find the run that achieved the highest ELBO
diags = stankmerdiag.collect()
process best_stankmer {
    input:
    file('*') from diags
    output:
    file('stan.best') into bestrun
    
    """
#!/usr/bin/env python
import sys, glob, re
max_d=1
max_r=-1
for file in glob.glob('stan.kmer.*.diag'):
    diagfile=open(file)
    m = re.search("kmer.(\\d+).diag",file)
    repid = m.group(1)
    for line in diagfile:
        if line[0] == '#': continue
        d=line.rstrip().split(",")
        if d[2] > max_d or max_d == 1:
            max_d = d[2]
            max_r = repid
sbfile = open("stan.best","w")
sbfile.write(max_r)
    """
}


// summarize the best run
runs = stankmerout.collect()
process summarize {
    publishDir 'out', mode: 'copy', overwrite: true

    input:
    file('unitigs.gfa') from gfa
    file('stan.kmer') from stankmer2
    file('stan.kmer.*.out') from runs
    file('stan.best') from bestrun

    output:
    file('strain_seqs.fa')
    """
    ${GDECONHOME}/stan_models/summarize.py stan.kmer stan.kmer.`cat stan.best`.out stan.kmer.summary ${params.relevance_threshold}
    ${GDECONHOME}/btools/consenseq.py unitigs.gfa stan.kmer.summary 0.51 strain_seqs.fa

    """
}
