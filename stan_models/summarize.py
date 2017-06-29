#!/usr/bin/env python
from __future__ import division
import sys
datafile = open(sys.argv[1])
postfile = open(sys.argv[2])
outfile = open(sys.argv[3], 'w')

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
samples = 0.0
for line in postfile:
    if line[0] == '#':
        continue
    if line[0] != '0':
        continue  # lp__ is always 0
    l = line.rstrip("\n").split(',')
    samples += 1.0
    # this is awful and hackish. parsing out node probabilities positionally instead of by looking at column headers
    i = len(l) - genomes*unitigs - 2*genomes
    for g in range(genomes):
        if not g in pprobs: pprobs[g] = {}
        for v in range(unitigs):
            if not v in pprobs[g]: pprobs[g][v] = 0.0
            pprobs[g][v] += float(l[i])
            i += 1
    # parse relevance in the same manner
    i = len(l) - 2*genomes
    for g in range(genomes):
        if not g in relevance: relevance[g] = 0.0
        relevance[g] += float(l[i])
        i+=1

for g in range(genomes):
    relevance[g] /= samples
    print "relevance " + str(relevance[g])

# write posterior mean
relevant=0
for v in range(unitigs):
    sep = ""
    for g in range(genomes):
        if relevance[g] > 1:
            continue # skip genome if not relevant
        relevant+=1
        outfile.write(sep)
        outfile.write( str(pprobs[g][v] / samples) )
        sep = "\t"
    outfile.write("\n")

print "found " + str(relevant/unitigs) + " genomes with " + str(unitigs) + " unitigs"

