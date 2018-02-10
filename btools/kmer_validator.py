#!/usr/bin/env python3
import gfapy
import sys
from Bio import SeqIO
from gfapy.sequence import rc

gfa = gfapy.Gfa.from_file(sys.argv[1])
klen = int(sys.argv[3])

refmers={}
graphmers={}

ref_kcount=0
for biorefseq in SeqIO.parse(sys.argv[2], 'fasta'):
    refseq = biorefseq.seq
    for i in range(len(refseq)):
        refmers[refseq[i:i+klen]]=1
        refmers[rc(refseq[i:i+klen])]=1
        ref_kcount+=1

print("Parsed "+str(ref_kcount)+" ref kmers, of which "+str(int(len(refmers)/2))+" are unique")


gfa_kcount=0
for seg in gfa.segments:
    for i in range(len(seg.sequence)):
        graphmers[seg.sequence[i:i+klen]]=1
        graphmers[rc(seg.sequence[i:i+klen])]=1
        gfa_kcount+=1
print("Parsed "+str(gfa_kcount)+" gfa kmers, of which "+str(int(len(graphmers)/2))+" are unique")


fn=0
tp=0
for kmer in refmers:
    if kmer in graphmers:
        tp+=1
        graphmers[kmer]=2 # mark as found
    else:
        fn+=1
fp=0
for kmer in graphmers:
    if graphmers[kmer]==1:
        fp+=1
tp=int(tp/2)
fp=int(fp/2)
fn=int(fn/2)

print("TP: "+str(tp)+" FP: "+str(fp)+" FN: "+str(fn))
