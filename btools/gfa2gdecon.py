#!/usr/bin/env python3
import gfapy
import sys
import pystan
import pickle
import numpy as np
import os


stan_source="unitig_reads_nmf.stan"

datatype = sys.argv[1]
gfa = gfapy.Gfa.from_file(sys.argv[2])
depthsfile = sys.argv[3]
linksfile = sys.argv[4]
strains = sys.argv[5]
min_g_size = sys.argv[6]
max_g_size = sys.argv[7]
repid = sys.argv[8]
known_covfile = sys.argv[9]
do_compile=False

sample_path = 'stan.'+repid+'.out'
diag_path = 'stan.'+repid+'.diag'

klen=0
if hasattr(gfa.header, 'kk') and gfa.header.kk is not None:
    klen = int(gfa.header.kk)

data={}

data['using_reads']=0
if datatype=='reads':
    data['using_reads']=1

# add basic info to data{}
data['G']=int(strains)
data['p_deadend']=0.00000000001
data['p_fork']=0.0000000001
data['genome_min']=int(min_g_size)
data['genome_max']=int(max_g_size)
data['V']=len(gfa.segments)
data['read_len']=150

dlabels = np.genfromtxt(depthsfile, delimiter=',', usecols=0, dtype=str)
if datatype=='reads':
    ddata = np.genfromtxt(depthsfile, delimiter=',', dtype=int)[:,1:]
    data['readcounts']=ddata
    data['kmercov']=np.zeros(ddata.shape)
else:
    ddata = np.genfromtxt(depthsfile, delimiter=',')[:,1:]
    data['kmercov']=ddata
    data['readcounts']=np.zeros(ddata.shape,dtype=int)

data['S']=ddata.shape[1]

if datatype=='reads':
    ldata = np.genfromtxt(linksfile, delimiter=',', usecols=(0,1), dtype=str)
    data['L']=ldata.shape[0]
    data['links']=np.ndarray(ldata.shape,dtype=int)
    lcdata = np.genfromtxt(linksfile, delimiter=',', dtype=int)[:,2:]
    data['linkcounts']=lcdata
else:
    data['L']=0
    data['links']=np.ndarray((0,2),dtype=int)
    data['linkcounts']=np.ndarray((0,data['S']),dtype=int)

cov_indices = np.genfromtxt(known_covfile,usecols=(0,1),dtype=int)
cov_vals = np.genfromtxt(known_covfile,usecols=4)
known_covs = np.zeros((data['S'],data['G']))
for i in range(cov_vals.shape[0]):
    known_covs[cov_indices[i,0]-1, cov_indices[i,1]-1]=cov_vals[i]

print(known_covs)
data['abund']=known_covs

segmap={}
for i in range(len(dlabels)):
    segmap[dlabels[i]]=i+1

if datatype=='reads':
    for i in range(ldata.shape[0]):
        if not ldata[i,0] in segmap: print("missing"+ldata[i,0])
        if not ldata[i,1] in segmap: print("missing"+ldata[i,1])
        data['links'][i,0]=segmap[ldata[i,0]]
        data['links'][i,1]=segmap[ldata[i,1]]

# unitig lengths
lenlist = np.ndarray((data['V']))
for seg in gfa.segments:
    lenlist[segmap[seg.name]-1] = len(seg.sequence)-klen
data['lengths']=lenlist

# parse and write out unitig adjacency structure
fromfirst = dict()
fromsecond = dict()
for edge in gfa.edges:
    fo, rfo = 1, -1
    if edge.from_orient == "-":
        fo, rfo = rfo, fo
    to, rto = 1, -1
    if edge.to_orient == "-":
        to, rto = rto, to

    # use segmap to map from segment names to array indices for the stan model
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

edgecounts = np.ndarray((data['V']))
for edge in gfa.edges:
    edgecounts[segmap[edge.from_segment.name]-1]+=1
    edgecounts[segmap[edge.to_segment.name]-1]+=1

data['coretigs']=(np.where(edgecounts==4)[0])+1
data['coretigcount']=len(data['coretigs'])
#print(data['coretigs'])

print(data)

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

fit = sm.vb(data=data, tol_rel_obj=0.001, output_samples=100, sample_file=sample_path, diagnostic_file=diag_path, algorithm="meanfield")
