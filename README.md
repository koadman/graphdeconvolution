# assembly graph deconvolution
A repository for work on deconvoluting assembly graphs containing strain mixtures using time-series abundance information.

### Quick start
Assuming you are starting in a directory containing a collection of paired-end Illumina read files, with names ending in the usual R?.fastq.gz, the software can be run as follows:

```
export BBMAP=/path/to/bbmap
export GDECONHOME=/path/to/this/repo
find `pwd` -maxdepth 1 -name "*R1.fastq.gz" | sort > read1_files.txt
find `pwd` -maxdepth 1 -name "*R2.fastq.gz" | sort > read2_files.txt
$GDECONHOME/btools/gdecon.nf --readlist1=read1_files.txt  --readlist2=read2_files.txt --r1suffix=R1.fastq.gz
```

If the workflow has completed successfully then the assembled strain genomes will appear in the directory `out/`.


### Dependencies and prerequisites

- Linux, kernel 2.6.32 or later
- nextflow, version 0.24 or later 
- python 3
- gfapy
- several others

### Major components of the workflow

The workflow first constructs a compacted de Bruijn graph using [bcalm](https://github.com/GATB/bcalm), then trims the dBg using [btrim](https://github.com/Malfoy/BTRIM). The abundance of each path (unitig) in each sample is then estimated via [tigops](https://github.com/GATB/tigops/).
Next, the abundance data and graph structure are given to a Bayesian graph deconvolution model.
The posterior output from the model is summarised and the number of strain genomes as well as their sequences is then recorded.

### Bayesian assembly graph deconvolution implemented in Stan

A model has been implemented in Stan code to carry out assembly graph deconvolution.
The input file for the model must contains the following information about the unitig graph:

- **V**: the number of unitigs
- **S**: the number of samples
- **depths**: a list of length VxS containing the depth of k-mer coverage for each unitig in each sample. Sample depths are given in succession for each unitig, e.g. if we denote the depth for unitig i in sample j as i.j, the depths are given as 1.1, 1.2, 1.3, ... , 2.1, 2.2, 2.3 and so on.
- **adj1count**: the number of unitigs with a single outgoing edge on an end
- **adj1source**: the source unitig ID for each unitig with a single outgoing edge on an end. Negative values indicate an outgoing edge coming from the end of the reverse-complement unitig sequence.
- **adj1dest**: the destination unitig ID corresponding to each of the source nodes listed in adj1source
- **adj2count**: the number of unitigs with two outgoing edges on an end
- **adj2source**: unitig IDs with two outgoing edges on an end
- **adj2dest1**: the first destination edge for the corresponding unitig given in adj2source
- **adj2dest2**: the second destination edge for the corresponding unitig given in adj2source
- **lengths**: the lengths of the V unitigs


### The MEGAHIT hack

We also created a hacked-up version of megahit, designed to take multiple samples. It uses the same command-line interface as standard megahit. Two additional output files are generated in addition to the usual megahit outputs. The first is `intermediate_contigs/k*.unitigs.fa`, which contains the unitig sequences. This can optionally be processed with the megahit toolkit contigs2fastg program to generate a graph viewable in BANDAGE. The second output is `intermediate_contigs/k*.unitig_depths.Rdata`. This file is as described above.

