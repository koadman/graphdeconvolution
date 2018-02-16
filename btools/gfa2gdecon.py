#!/usr/bin/env python3
import gfapy
import sys
import pystan
import pickle
import numpy as np
import os

stan_source="genotypes3_unitig_graph.stan"

gfa = gfapy.Gfa.from_file(sys.argv[1])
depthsfile = sys.argv[2]
strains = sys.argv[3]
min_g_size = sys.argv[4]
max_g_size = sys.argv[5]
repid = sys.argv[6]
do_compile=False

sample_path = 'stan.'+repid+'.out'
diag_path = 'stan.'+repid+'.diag'

klen = int(gfa.header.kk)

data={}

# add basic info to data{}
data['G']=int(strains)
data['p_deadend']=0.0000001
data['p_fork']=0.00001
data['genome_min']=int(min_g_size)
data['genome_max']=int(max_g_size)
data['V']=len(gfa.segments)

dlabels = np.genfromtxt(depthsfile, delimiter=',', usecols=0, dtype=str)
ddata = np.genfromtxt(depthsfile, delimiter=',')[:,1:]
data['depths']=ddata
data['S']=ddata.shape[1]

segmap={}
for i in range(len(dlabels)):
    segmap[dlabels[i]]=i+1

# unitig lengths
lenlist = []
for seg in gfa.segments:
    lenlist.append(len(seg.sequence)-klen)
data['lengths']=lenlist

# parse and write out unitig adjacency structure
fromfirst = dict()
fromsecond = dict()
for edge in gfa.edges:
    fo = 1
    rfo = -1
    if edge.from_orient == "-":
        fo = -1
        rfo = 1
    to = 1
    rto = -1
    if edge.from_orient == "-":
        to = -1
        rto = 1

    # assume integer unitig IDs indexed from 0 -- and convert these to 1-based
    fname = fo*segmap[edge.from_segment.name]
    tname = to*segmap[edge.to_segment.name]
    if fname not in fromfirst:
        fromfirst[fname] = tname
    else:
        fromsecond[fname] = tname
    rfname = rfo*segmap[edge.from_segment.name]
    rtname = rto*segmap[edge.to_segment.name]
    if rtname not in fromfirst:
        fromfirst[rtname] = rfname
    else:
        fromsecond[rtname] = rfname


data['adj1count']=len(fromfirst)-len(fromsecond)
data['adj2count']=len(fromsecond)

data['adj2source']=[]
data['adj2dest1']=[]
data['adj2dest2']=[]
data['adj1source']=[]
data['adj1dest']=[]

for seg in fromfirst:
    if seg in fromsecond:
        data['adj2source'].append(seg)
        data['adj2dest1'].append(fromfirst[seg])
        data['adj2dest2'].append(fromsecond[seg])
    else:
        data['adj1source'].append(seg)
        data['adj1dest'].append(fromfirst[seg])


# Path to the compiled stan binary
my_path = os.path.split(os.path.realpath(__file__))[0]
binary = os.path.join(my_path, '..', 'stan_models', stan_source + '.pkl')
source_md5 = os.path.join(my_path, '..', 'stan_models', stan_source + '.md5')
source_file  = os.path.join(my_path, '..', 'stan_models', stan_source)

# Compile a stan binary if it does not exist or compile option is set
if do_compile or not os.path.isfile(binary):
    sm = pystan.StanModel(file=source_file)
    with open(binary, 'wb') as f:
        pickle.dump(sm, f)
else:
    sm = pickle.load(open(binary, 'rb'))

fit = sm.vb(data=data, tol_rel_obj=0.008, sample_file=sample_path, diagnostic_file=diag_path, algorithm="meanfield")
