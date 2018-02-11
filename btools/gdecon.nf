#!/usr/bin/env nextflow
/**
 * Usage: gdecon.nf --readlist1=/path/to/pe_readfile_list.txt  --readlist2=/path/to/pe_readfile_list.txt --r1suffix=.r1.fastq.gz
 * (c) 2017 Aaron Darling
 */

// these default settings can be overridden on the command-line
params.mingsize=100000
params.maxgsize=600000
params.k = 51
params.maxstrains = 25
params.relevance_threshold = 0
params.maxtiplen = 100
params.asmcores = 8
params.covk = 51
params.output="out"

rfiles1c = Channel.from(file(params.readlist1)).splitText() { it.trim() }
rfiles2c = Channel.from(file(params.readlist2)).splitText() { it.trim() }
rfiles_cleaning = rfiles1c.merge(rfiles2c) {o, e -> [o, e]}

process clean {
    input:
    set r1,r2 from rfiles_cleaning
    output:
    file('*.clean.fq') into cleanreads

    """
    name=`basename ${r1} ${params.r1suffix}`
    ${BBMAP}/bbduk.sh -Xmx500m in=${r1} in2=${r2} out=\$name.clean.fq outs=\$name.single.fq ref=${BBMAP}/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 qtrim=r trimq=20 tpe tbo
    """
}

// constructs the unitig graph
allclean = cleanreads.collect()
process unitig {
    input:
    file('*') from allclean
    output:
    file('rlist.contigs.fa') into unitigs
    file('rlist.contigs.fa') into unitigs2

    """
    ls -1 *.clean.fq > rlist
    ${GDECONHOME}/external/gatb/minia -in rlist -kmer-size 51 -no-bulge-removal -ec-rctc-cutoff 40 -tip-rctc-cutoff 10 -abundance-min 3 -abundance-min-threshold 3
#    -tip-len-topo-kmult 3 -tip-len-rctc-kmult 3 -tip-rctc-cutoff 15 -bulge-len-kmult 3 -bulge-altpath-covmult 20 -bulge-altpath-kadd 1000 -bulge-len-kadd 52 -ec-rctc-cutoff 30 -ec-len-kmult 1.5
# -in rlist -kmer-size 51 -tip-rctc-cutoff 10 -bulge-len-kmult 2 -bulge-altpath-covmult 50 -bulge-len-kadd 52 -ec-rctc-cutoff 30 -ec-len-kmult 2
# -in rlist -kmer-size ${params.k} -tip-rctc-cutoff 10 -bulge-len-kmult 2 -bulge-altpath-covmult 10 -bulge-len-kadd 50
    """
}

rfiles1 = Channel.from(file(params.readlist1)).splitText() { it.trim() }
rfiles2 = Channel.from(file(params.readlist2)).splitText() { it.trim() }
rfiles = rfiles1.merge(rfiles2) {o, e -> [o, e]}.spread(unitigs)

// estimates depth of coverage on unitigs
process unitig_cov {
    input:
    set reads,reads2,file(unitigs) from rfiles
    output:
    file('*.counts') into covcounts

    """
    name=`basename $reads ${params.r1suffix}`
    cat $reads $reads2 > bothreads
    ${GDECONHOME}/external/gatb/tigops coverage -out \$name.counts -kmer-size ${params.covk} -reads bothreads -tigs $unitigs -name \$name
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
        m=re.search('_cov_([^_]+?)_', line)
        if m == None: print "line " + line + " in " + filename + " contains unexpected formatting"
        if not ti in tigs: tigs[ti]=str(ti)
        tigs[ti] += "," + str(float(m.group(1)))
        ti+=1
for i in range(ti):
    covout.write(tigs[i]+"\\n")

"""
}

// convert the bcalm unitigs to GFA format
process convert_gfa {
    input:
    file('unitigs.fa') from unitigs2

    output:
    file('unitigs.gfa') into gfa
    file('unitigs.gfa') into gfa2

"""
    ${GDECONHOME}/external/gatb/convertToGFA.py unitigs.fa unitigs.gfa ${params.k}
    perl -p -i -e "s/\\s+k:i:/\\tkk:i:/g" unitigs.gfa
#    perl -p -i -e "s/KM:f:/km:f:/g" unitigs.gfa
"""
}

// run the Bayesian NMF in stan
process stankmer_nmf {
    input:
    each repid from 1,2,3,4
    file('unitigs.gfa') from gfa
    file('coverage.csv') from coverage
    output:
    file("stan.*.diag") into stankmerdiag
    file("stan.*.out") into stankmerout

    """
    ${GDECONHOME}/btools/gfa2gdecon.py unitigs.gfa coverage.csv ${params.maxstrains} ${params.mingsize} ${params.maxgsize} ${repid}
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
max_d=1.0
max_r=-1
for file in glob.glob('stan.*.out.diag'):
    diagfile=open(file)
    m = re.search("(\\d+).out.diag",file)
    repid = m.group(1)
    for line in diagfile:
        if line[0] == '#': continue
        d=line.rstrip().split(",")
        if float(d[2]) > max_d or max_d == 1:
            max_d = float(d[2])
            max_r = repid
sbfile = open("stan.best","w")
sbfile.write(max_r)
    """
}


// select the best run of variational inference
runs = stankmerout.collect()
process best_stankmer {
    input:
    file('stan.*.out') from runs
    file('stan.best') from bestrun

    output:
    file('stan.kmer.out') into beststankmer
    file('stan.kmer.out') into beststankmer2

    """
    cp stan.`cat stan.best`.out stan.kmer.out
    """
}

// summarize the ARD weights
process get_ard_weights {
    publishDir params.output, mode: 'copy', overwrite: true

    input:
    file('unitigs.gfa') from gfa
    file('stan.kmer.out') from beststankmer

    output:
    file('ard_weights.txt') into relweights

    when:
    params.relevance_threshold == 0

    script:
    """
     ${GDECONHOME}/stan_models/summarize.py unitigs.gfa stan.kmer.out stan.kmer.summary ${params.maxstrains} 0 | sort -n > ard_weights.txt
    """
}

relweights.subscribe {
    println("\nStrain count evaluation complete. Please inspect the file ${params.output}/ard_weights.txt and identify a numeric value that separates the weights into two natural groups.\nThen re-run gdecon.nf with the following command, specifying the chosen threshold with --relevance_threshold:\n")
    println(workflow.commandLine + " --relevance_threshold=<number> ")
}

// summarize the best run
process summarize {
    publishDir params.output, mode: 'copy', overwrite: true

    input:
    file('unitigs.gfa') from gfa2
    file('stan.kmer.out') from beststankmer2

    output:
    file('strain_seqs.fa') into strainseqs

    when:
    params.relevance_threshold != 0

    script:
    """
    ${GDECONHOME}/stan_models/summarize.py unitigs.gfa stan.kmer.out stan.kmer.summary ${params.maxstrains} ${params.relevance_threshold}
    ${GDECONHOME}/btools/consenseq.py unitigs.gfa stan.kmer.summary 0.5 strain_seqs.fa
    """
}

strainseqs.subscribe {
    println("\nWorkflow complete. Assembled strain sequences have been stored in ${params.output}/strain_seqs.fa")
}
