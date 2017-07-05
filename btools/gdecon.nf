#!/usr/bin/env nextflow
/**
 * Usage: gdecon.nf --readlist1=/path/to/pe_readfile_list.txt  --readlist2=/path/to/pe_readfile_list.txt --r1suffix=.r1.fastq.gz
 * (c) 2017 Aaron Darling
 */

readlist1 = Channel.from(file(params.readlist1))
readlist2 = Channel.from(file(params.readlist2))

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
    ${GDECONHOME}/external/gatb/bcalm -in bothreads.txt -kmer-size 61 -abundance-min 2 -nb-cores 8
    ${GDECONHOME}/external/gatb/btrim bothreads.unitigs.fa 61 100 8
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
    ${GDECONHOME}/external/gatb/tigops coverage -out \$name.counts -kmer-size 31 -reads bothreads -tigs unitigs -name \$name
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
sep=""
for filename in listing:
    ti = 0
    f = open(filename)
    for line in f:
        if line[0] != '>':
            continue
        m=re.search('mean_(\\d*\\.*\\d*)_', line)
# _mean_17.8_med
        if not ti in tigs: tigs[ti]=""
        tigs[ti] += sep + m.group(1)
        ti+=1
    sep=","
for i in range(ti):
    covout.write(tigs[i]+"\\n")

"""
}

// run the Bayesian NMF in stan
process stankmer_nmf {
    publishDir 'out', mode: 'copy', overwrite: true

    input:
    file('coverage.csv') from coverage
    file('unitigs.fa') from unitigs2

    output:
    file('strain_seqs.fa')

    """
    ${GDECONHOME}/external/gatb/convertToGFA.py unitigs.fa unitigs.gfa 51
    perl -p -i -e "s/\\s+k:i:/\\tkk:i:/g" unitigs.gfa
    perl -p -i -e "s/KM:f:/km:f:/g" unitigs.gfa
    ${GDECONHOME}/btools/gfa2gdecon.py unitigs.gfa coverage.csv 3 1500 > stan.kmer
    ${GDECONHOME}/stan_models/genotypes3_unitig_graph variational output_samples=100 tol_rel_obj=0.008 iter=10000 data file=stan.kmer output file=stan.kmer.out diagnostic_file=stan.kmer.diag
    ${GDECONHOME}/stan_models/summarize.py stan.kmer stan.kmer.out stan.kmer.summary 2
    ${GDECONHOME}/btools/consenseq.py unitigs.gfa stan.kmer.summary 0.5 strain_seqs.fa
    """
}
