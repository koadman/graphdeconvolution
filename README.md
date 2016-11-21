# assembly graph deconvolution
A repository for work on deconvoluting assembly graphs containing strain mixtures using time-series abundance information.

The code structure consists of two major parts:
1. A hacked-up version of MEGAHIT which produces unitigs and associated depth of coverage information.
2. A deconvolution algorithm (or many of these as we test out different approaches)

### The MEGAHIT hack

The hacked-up version of megahit uses the same command-line interface as standard megahit. Two additional output files are generated in addition to the usual megahit outputs. The first is `intermediate_contigs/k*.unitigs.fa`, which contains the unitig sequences. This can optionally be processed with the megahit toolkit contigs2fastg program to generate a graph viewable in BANDAGE. The second output is `intermediate_contigs/k*.unitig_depths.Rdata`. This file contains the following information about the unitig graph:

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


### Step-by-step instructions

```
cd megahit
make
cd ../test_data
../megahit/megahit -1 `cat testr1A.csv` -2 `cat r2A.csv` -o asmtest

# specify the number of genomes (factors) as 4 and the maximum depth of coverage
echo "G<-4" > unitigs.Rdata
echo "maxdepth<-200" >> unitigs.Rdata
cat asmtest/intermediate_contigs/k99.unitig_depths.Rdata >> unitigs.Rdata

# run the variational inference
../stan_models/genotypes3_unitig_graph  variational output_samples=100 tol_rel_obj=0.001 iter=10000 data file=unitigs.Rdata output file=factor_16S.4.out diagnostic_file=factor_16S.4.diag

# the following command summarizes the posterior probability of each unitig being in each of the (four) genomes
../stan_models/summarize.py unitigs.Rdata factor_16S.4.out factor_16S.4.summary

# the next command parses the unitig fasta and the posterior summary and creates one FastG file per genome with the base name 'fourbins'
# it uses a posterior threshold of 0.5 for contig inclusion
../megahit/megahit_toolkit binUnitigs 99 asmtest/intermediate_contigs/k99.unitigs.fa factor_16S.4.summary fourbins 0.5
```

