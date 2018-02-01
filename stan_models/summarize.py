#!/usr/bin/env python
from __future__ import division
import sys
import numpy as np
if len(sys.argv) != 5:
    print "Usage: summarize.py <stan data> <stan posterior> <FastA output> <relevance threshold>"
    sys.exit(-1)
datafile = open(sys.argv[1])
postfile = open(sys.argv[2])
outfile = open(sys.argv[3], 'w')
relevance_threshold = float(sys.argv[4])


genomes = 0
unitigs = 0
for line in datafile:
    line = line.rstrip("\n")
    if line.startswith("G<-"):
        genomes = int(line[3:])
    if line.startswith("V<-"):
        unitigs = int(line[3:])


pprobs = {}
relevance = {}

popofile = np.genfromtxt(sys.argv[2],comments='#',skip_header=28,delimiter=',',names=True,deletechars="""~!@#$%^&*()=+~\|]}[{';: /?>,<""")

for u in range(unitigs):    
    for g in range(genomes):
        if not g in pprobs: pprobs[g] = {}
        genotxt=popofile['geno.'+str(u+1)+'.'+str(g+1)]        
        pprobs[g][u]=np.sum(genotxt)/len(genotxt)

for r in range(genomes):
    rtxt=popofile['relevance.'+str(r+1)]
    relevance[r]=np.sum(rtxt)/len(rtxt)
    print("relevance "+str(relevance[r]))

# write posterior mean
relevant=0
for v in range(unitigs):
    sep = ""
    for g in range(genomes):
        if relevance[g] > relevance_threshold:
            continue # skip genome if not relevant
        relevant+=1
        outfile.write(sep)
        outfile.write( str(pprobs[g][v]) )
        sep = "\t"
    outfile.write("\n")

print "found " + str(relevant/unitigs) + " genomes with " + str(unitigs) + " unitigs"

