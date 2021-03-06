#!/usr/bin/env nextflow
/**
 * Usage: gdecon.nf --readlist1=/path/to/pe_readfile_list.txt  --readlist2=/path/to/pe_readfile_list.txt --r1suffix=.r1.fastq.gz
 * (c) 2017 Aaron Darling
 */

// these default settings can be overridden on the command-line
params.mingsize=100000
params.maxgsize=600000
params.k = 51
params.maxstrains = 4
params.relevance_threshold = 0
params.maxtiplen = 100
params.asmcores = 8
params.covk = params.k
params.output="out"

known_coverages = Channel.from(file(params.coverages))

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
    ${GDECONHOME}/external/gatb/minia -in rlist -kmer-size ${params.k} -no-bulge-removal -ec-rctc-cutoff 40 -tip-rctc-cutoff 10 -abundance-min 3 -abundance-min-threshold 3
    rm -f *glue*
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
    file('coverage.csv') into coverage2

"""
#!/usr/bin/env python
import glob
import re
listing = glob.glob('*.counts')
covout = open('coverage.csv', 'w')
tigs = {}
#for filename in listing:
#    f = open(filename)
for i in range(16):
    f = open("seed#1234-+-alpha#0.1-+-xfold#30.wgs."+str(i+1)+"..counts")
    for line in f:
        if line[0] != '>':
            continue
        m=re.search('_cov_([^_]+?)_', line)
        if m == None: print "line " + line + " in " + filename + " contains unexpected formatting"
        ti = line.split()[0][1:]
        if not ti in tigs:
          tigs[ti]=ti
        tigs[ti] += "," + str(float(m.group(1)))
for ti in tigs:
    covout.write(tigs[ti]+"\\n")

"""
}

// convert the bcalm unitigs to GFA format
process convert_gfa {
  publishDir params.output, mode: 'copy', overwrite: true
    input:
    file('unitigs.fa') from unitigs2

    output:
    file('unitigs.gfa') into gfa
    file('unitigs.gfa') into gfa2
    file('unitigs.gfa') into gfa3

"""
    ${GDECONHOME}/external/gatb/convertToGFA.py unitigs.fa unitigs.gfa ${params.k}
    perl -p -i -e "s/\\s+k:i:/\\tkk:i:/g" unitigs.gfa
"""
}

// convert the bcalm unitigs to GFA format
process vgify_asm {
    publishDir params.output, mode: 'copy', overwrite: true
    input:
    file('unitigs.gfa') from gfa3

    output:
    file('split.vg') into vg_asm
    file('split.gfa') into vg_gfa
    file('split.gfa') into vg_gfa2

"""
    vg view -vF unitigs.gfa > unitigs.vg
    vg mod -U 10 unitigs.vg > unitigs.norm.vg
    vg mod -X 1024 unitigs.norm.vg > split.vg
    vg view -g split.vg > split.gfa
    vg view -j split.vg > split.json
    vg index split.vg -x split.vg.xg -g split.vg.gcsa -k 9
"""
}

vgrfiles1 = Channel.from(file(params.readlist1)).splitText() { it.trim() }
vgrfiles2 = Channel.from(file(params.readlist2)).splitText() { it.trim() }
vgrfiles = vgrfiles1.merge(vgrfiles2) {o, e -> [o, e]}.spread(vg_asm)

process map_reads {
    input:
    set reads,reads2,file('split.vg') from vgrfiles
    output:
    file('*.counts') into vgcovcounts
    file('*.links') into vglinkcounts

"""
    vg index split.vg -x split.vg.xg -g split.vg.gcsa -k 9
    vg map -k 9 -d split.vg -x split.vg.xg -g split.vg.gcsa -f $reads > aligned1.gam
    vg map -d split.vg -x split.vg.xg -g split.vg.gcsa -f $reads2 > aligned2.gam
    vg view -aj aligned1.gam >> aligned.json
    vg view -aj aligned2.gam >> aligned.json
    vg view -j split.vg > split.json
    name=`basename $reads ${params.r1suffix}`
     ${GDECONHOME}/btools/node_cov.py split.json aligned.json \$name.counts \$name.links
"""
}

vgcovfefe = vgcovcounts.collect()
process vg_extract_cov {
    input:
    file('*') from vgcovfefe
    output:
    file('coverage.csv') into vgcoverage
    file('coverage.csv') into vgcoverage2

"""
#!/usr/bin/env python
import glob
import re
listing = glob.glob('*.counts')
covout = open('coverage.csv', 'w')
tigs = {}
#for filename in listing:
#    f = open(filename)
for i in range(16):
    f = open("seed#1234-+-alpha#0.1-+-xfold#30.wgs."+str(i+1)+"..counts")
    for line in f:
        (ti,cov) = line.rstrip().split(',')
        if not ti in tigs:
          tigs[ti]=ti
        tigs[ti] += "," + cov
for ti in tigs:
    covout.write(tigs[ti]+"\\n")
"""
}

vglinky = vglinkcounts.collect()
process vg_extract_links {
    input:
    file('*') from vglinky
    output:
    file('links.csv') into vglinks

"""
#!/usr/bin/env python
import glob
import re
listing = glob.glob('*.links')
tigs = {}
f_i=-1
#for filename in listing:
#    f = open(filename)
for i in range(16):
    f = open("seed#1234-+-alpha#0.1-+-xfold#30.wgs."+str(i+1)+"..links")
    f_i+=1
    for line in f:
        (t1,t2,links) = line.rstrip().split('\t')
        if not (t1,t2) in tigs:
          tigs[(t1,t2)]=['0']*len(listing)
        tigs[(t1,t2)][f_i] = links
covout = open('links.csv', 'w')
for ti in tigs:
    covout.write(",".join(ti)+','+",".join(tigs[ti])+"\\n")
"""
}


// run the Bayesian NMF in stan
process stankmer_nmf_reads {
    input:
    each repid from 1,2,3,4
    file('unitigs_vg.gfa') from vg_gfa
    file('read_coverage.csv') from vgcoverage
    file('links.csv') from vglinks
    file('known_coverage.tsv') from known_coverages

    file('unitigs.gfa') from gfa
    file('coverage.csv') from coverage
    output:
    file("stan.*.diag") into stankmerdiag
    file("stan.*.out") into stankmerout

    """
    ${GDECONHOME}/btools/gfa2gdecon.py reads unitigs_vg.gfa read_coverage.csv links.csv ${params.maxstrains} ${params.mingsize} ${params.maxgsize} ${repid} known_coverage.tsv
#    ${GDECONHOME}/btools/gfa2gdecon.py kmercov unitigs.gfa coverage.csv links.csv ${params.maxstrains} ${params.mingsize} ${params.maxgsize} ${repid} known_coverage.tsv
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
for file in glob.glob('stan.*.diag'):
    diagfile=open(file)
    m = re.search("(\\d+).diag",file)
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
    file('stan.kmer.out') from beststankmer
    file('coverage.csv') from vgcoverage2

    output:
    file('ard_weights.txt') into relweights

    when:
    params.relevance_threshold == 0

    script:
    """
     ${GDECONHOME}/stan_models/summarize.py coverage.csv stan.kmer.out stan.kmer.summary ${params.maxstrains} 0 | sort -n > ard_weights.txt
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
    file('unitigs_vg.gfa') from vg_gfa2
    file('stan.kmer.out') from beststankmer2
    file('coverage.csv') from coverage2
    file('coverage_vg.csv') from vgcoverage2

    output:
    file('assembled.*') into strainseqs

    when:
    params.relevance_threshold != 0

    script:
    """
#    ${GDECONHOME}/stan_models/summarize.py coverage.csv stan.kmer.out stan.kmer.summary ${params.maxstrains} ${params.relevance_threshold}
#    ${GDECONHOME}/btools/consenseq.py unitigs.gfa coverage.csv stan.kmer.summary 0.5 assembled
    ${GDECONHOME}/stan_models/summarize.py coverage_vg.csv stan.kmer.out stan.kmer.summary ${params.maxstrains} ${params.relevance_threshold}
    ${GDECONHOME}/btools/consenseq.py unitigs_vg.gfa coverage_vg.csv stan.kmer.summary 0.5 assembled
    """
}

strainseqs.subscribe {
    println("\nWorkflow complete. Assembled strain sequences have been stored in ${params.output}/assembled.*.fasta")
}
