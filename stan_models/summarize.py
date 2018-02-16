#!/usr/bin/env python3
from __future__ import division
import sys
import gfapy
import numpy as np
if len(sys.argv) != 6:
    print("Usage: summarize.py <assembly GFA> <stan posterior> <posterior summary output> <relevance threshold>")
    sys.exit(-1)

gfa = gfapy.Gfa.from_file(sys.argv[1])
postfile = sys.argv[2]
outfile = open(sys.argv[3], 'w')
genomes = int(sys.argv[4])
relevance_threshold = float(sys.argv[5])

unitigs = len(gfa.segments)

pprobs = {}
relevance = {}
glen = {}

with open(postfile) as f:
    lines = (line for line in f if not line.startswith('#'))
    popofile = np.genfromtxt(lines,comments='#',delimiter=',',names=True,deletechars="""~!@#$%^&*()=+~\|]}[{';: /?>,<""")

for u in range(unitigs):
    for g in range(genomes):
        if not g in pprobs: pprobs[g] = {}
        genotxt=popofile['geno.'+str(u+1)+'.'+str(g+1)]
        pprobs[g][u]=np.sum(genotxt)/len(genotxt)

for r in range(genomes):
    rtxt=popofile['relevance.'+str(r+1)]
    relevance[r]=np.sum(rtxt)/len(rtxt)
    print("relevance "+str(relevance[r]))

glentxt=popofile['genome_size']
glen[1]=np.sum(glentxt)/len(glentxt)
print("length "+str(glen[1]))

# write posterior mean
relevant=0
for g in range(genomes):
    if relevance[g] > relevance_threshold:
        continue # skip genome if not relevant
    outfile.write("\tgenome_"+str(relevant))
    relevant+=1
if relevant==0: outfile.write("\t")
outfile.write("\n")

for v in range(unitigs):
    outfile.write(gfa.segment_names[v])
    for g in range(genomes):
        if relevance[g] > relevance_threshold:
            continue # skip genome if not relevant
        outfile.write( "\t"+str(pprobs[g][v]) )
    outfile.write("\n")

print("found " + str(relevant/unitigs) + " genomes with " + str(unitigs) + " unitigs")
