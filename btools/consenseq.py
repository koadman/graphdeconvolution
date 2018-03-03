#!/usr/bin/env python3
import gfapy
import sys
import numpy as np

gfa = gfapy.Gfa.from_file(sys.argv[1])
covfile = sys.argv[2]
summaryfile = sys.argv[3]
posterior_threshold = float(sys.argv[4])
klen=1
if hasattr(gfa.header, 'kk') and gfa.header.kk is not None:
    klen = int(gfa.header.kk)

segmap={}
invmap={}
dlabels = np.genfromtxt(summaryfile, usecols=0, skip_header=1, dtype=str)
for i in range(len(dlabels)):
    segmap[dlabels[i]]=i
    invmap[i]=dlabels[i]

print("Parsing posterior summary")
posts = np.genfromtxt(summaryfile, skip_header=1)[:,1:]
visited = np.zeros(posts.shape)
print("Parsed "+str(visited.shape[1])+" strain posteriors")

strainseqs = ["" for s in range(posts.shape[1])]
strainpaths = ["" for s in range(posts.shape[1])]
segs = gfa.segments
seqmap={}
for seg in gfa.segments:
    seqmap[segmap[seg.name]]=seg.sequence.rstrip()

edgemap=[[] for i in range(len(dlabels))]
for e in gfa.edges:
    edgemap[segmap[e.from_segment.name]].append(e)
    edgemap[segmap[e.to_segment.name]].append(e)

for s in range(posts.shape[1]):
    pathcount = 0
    for n in range(posts.shape[0]):
        if visited[n,s] == 1:
            continue
        if posts[n,s] < posterior_threshold:
            continue
        # this node exists in the strain and has not yet been visited
        # start a graph traversal in both directions
        strainpath = "+"+str(invmap[n])
        cur_seq = seqmap[n]
        cur_seg = n
        cur_orient = "+"
        visited[n,s] = 1
        for i in range(2):
            while True:
                # find possible successors of the current segment
                successors = dict()
                if edgemap[cur_seg]==0:
                    print("isolated node "+invmap[cur_seg])
                for e in edgemap[cur_seg]:
                    if e.from_segment.name == invmap[cur_seg] and e.from_orient == cur_orient and posts[segmap[e.to_segment.name],s] >= posterior_threshold:
                        successors[segmap[e.to_segment.name]] = e.to_orient
                    if e.to_segment.name == invmap[cur_seg] and e.to_orient != cur_orient and posts[segmap[e.from_segment.name],s] >= posterior_threshold:
                        successors[segmap[e.from_segment.name]] = "-" if e.from_orient == "+" else "+"

                if len(successors) > 1:
                    print("ambiguous successor for node " + invmap[cur_seg] + " strain " + str(s))
                ambig = False
                for suc in successors:
                    if visited[suc,s] == 1:
                        ambig = True    # if the node was already visited, and the path was not extended, then it was likely ambiguous (or a repeat?)
                if len(successors) != 1 or ambig:
                    break   ## the successor node is either ambiguous or nonexistent
                for suc in successors:
                    if visited[suc,s] == 1:
                        print("Error already visited node " + invmap[suc] + " in strain " + str(s))
                    cur_orient = successors[suc]
                    cur_seg = suc
                    visited[suc,s] = 1
                    if i == 0:
                        strainpath += "," + cur_orient + invmap[suc]
                    else:
                        strainpath = cur_orient + invmap[suc] + "," + strainpath
                    if(cur_orient == "+"):
                        cur_seq += seqmap[suc][klen-1:]
                    else:
                        cur_seq += gfapy.sequence.rc(seqmap[suc])[klen-1:]

            # now try to extend in the other direction from the original node
            cur_seg = n
            cur_orient = "-"
            cur_seq = gfapy.sequence.rc(cur_seq)

        strainseqs[s] += ">strain_" + str(s) + "_path_" + str(pathcount) + "\n" + cur_seq + "\n"
        strainpaths[s] += ">strain_" + str(s) + "_path_" + str(pathcount) + "\n" + strainpath + "\n"
        pathcount += 1

i=0
for s in strainseqs:
    outseqs = open(sys.argv[5] + ".strain_"+str(i)+".fasta","w")
    outseqs.write(s)
    outseqs.close()
    i+=1

for s in strainpaths:
    print(s.rstrip())
