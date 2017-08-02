#!/usr/bin/env python3
import gfapy
import sys

print("Parsing GFA")
gfa = gfapy.Gfa.from_file(sys.argv[1])
print("Done parsing GFA")
summary = open(sys.argv[2])
posterior_threshold = float(sys.argv[3])
klen = gfa.header.kk

print("making segmap")
segmap = dict()
segnames = gfa.segment_names
for seg in range(len(segnames)):
    segmap[int(segnames[seg])]=seg
print("done making segmap")

print("Parsing posterior summary")
posts = list()
visited = list()
for line in summary:
    dd = line.rstrip().split()
    for i in range(len(dd)):
        if len(posts) == 0:
            posts = [[] for j in range(len(dd))]
            visited = [[] for j in range(len(dd))]
        posts[i].append(float(dd[i]))
        visited[i].append(0)

print("Parsed "+str(len(visited))+" strain posteriors")

strainseqs = ["" for s in range(len(posts))]
strainpaths = ["" for s in range(len(posts))]
segs = gfa.segments
edgemap=dict()
for e in gfa.edges:
    fname = int(e.from_segment.name)
    tname = int(e.to_segment.name)
    if not fname in edgemap:
        edgemap[fname]=[]
    if not tname in edgemap:
        edgemap[tname]=[]
    edgemap[fname].append(e)
    edgemap[tname].append(e)

for s in range(len(posts)):
    pathcount = 0
    for n in range(len(posts[s])):
        if visited[s][n] == 1:
            continue
        if posts[s][n] < posterior_threshold:
            continue
        # this node exists in the strain and has not yet been visited
        # start a graph traversal in both directions
        strainpath = "+"+str(n)
        cur_seq = segs[segmap[n]].sequence.rstrip()
        cur_seg = n
        cur_orient = "+"
        visited[s][n] = 1
        for i in range(2):
            while True:
                # find possible successors of the current segment
                successors = dict()
                if not cur_seg in edgemap:
                    edgemap[cur_seg]=[]
                for e in edgemap[cur_seg]:
                    if e.from_segment.name == str(cur_seg) and e.from_orient == cur_orient and posts[s][int(e.to_segment.name)] >= posterior_threshold:
                        successors[int(e.to_segment.name)] = e.to_orient
                    if e.to_segment.name == str(cur_seg) and e.to_orient != cur_orient and posts[s][int(e.from_segment.name)] >= posterior_threshold:
                        successors[int(e.from_segment.name)] = "-" if e.from_orient == "+" else "+"

                if len(successors) > 1:
                    print("ambiguous successor for node " + str(cur_seg) + " strain " + str(s))
                ambig = False
                for suc in successors:
                    if visited[s][suc] == 1:
                        ambig = True    # if the node was already visited, and the path was not extended, then it was likely ambiguous (or a repeat?)
                if len(successors) != 1 or ambig:
                    break   ## the successor node is either ambiguous or nonexistent
                for suc in successors:
                    if visited[s][suc] == 1:
                        print("Error already visited node " + str(suc) + " in strain " + str(s))
                    cur_orient = successors[suc]
                    cur_seg = suc
                    visited[s][suc] = 1
                    if i == 0:
                        strainpath += "," + cur_orient + str(suc)
                    else:
                        strainpath = cur_orient + str(suc) + "," + strainpath
                    if(cur_orient == "+"):
                        cur_seq += segs[segmap[suc]].sequence[klen-1:]
                    else:
                        cur_seq += gfapy.sequence.rc(segs[segmap[suc]].sequence)[klen-1:]

            # now try to extend in the other direction from the original node
            cur_seg = n
            cur_orient = "-"
            cur_seq = gfapy.sequence.rc(cur_seq)

        strainseqs[s] += ">strain_" + str(s) + "_path_" + str(pathcount) + "\n" + cur_seq + "\n"
        strainpaths[s] += ">strain_" + str(s) + "_path_" + str(pathcount) + "\n" + strainpath + "\n"
        pathcount += 1

outseqs = open(sys.argv[4],"w")
for s in strainseqs:
    outseqs.write(s)

for s in strainpaths:
    print(s.rstrip())
