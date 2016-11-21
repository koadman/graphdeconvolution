#!/usr/bin/env python
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

print "found " + str(genomes) + " genomes with " + str(unitigs) + " unitigs"

pprobs = {}

samples = 0.0
for line in postfile:
    if line[0] == '#':
        continue
    if line[0] != '0':
        continue  # lp__ is always 0
    l = line.rstrip("\n").split(',')
    samples += 1.0
    i = len(l) - genomes*unitigs - genomes
    for g in range(genomes):
        if not g in pprobs: pprobs[g] = {}
        for v in range(unitigs):
            if not v in pprobs[g]: pprobs[g][v] = 0.0
            pprobs[g][v] += float(l[i])
            i += 1


# write posterior mean
for v in range(unitigs):
    sep = ""
    for g in range(genomes):
        outfile.write(sep)
        outfile.write( str(pprobs[g][v] / samples) )
        sep = "\t"
    outfile.write("\n")

