#!/usr/bin/env python3
import gfapy
import sys

gfa = gfapy.Gfa.from_file(sys.argv[1])
summary = open(sys.argv[2])
posterior_threshold = float(sys.argv[3])
klen = gfa.header.kk

segmap = dict()
for seg in range(len(gfa.segments)):
    segmap[int(gfa.segments[seg].name)]=seg

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


strainseqs = ["" for s in range(len(posts))]
strainpaths = ["" for s in range(len(posts))]
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
        cur_seq = gfa.segments[segmap[n]].sequence.rstrip()
        cur_seg = n
        cur_orient = "+"
        visited[s][n] = 1
        for i in range(2):
            while True:
                # find possible successors of the current segment
                successors = dict()
                for e in gfa.edges:
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
                        cur_seq += gfa.segments[segmap[suc]].sequence[klen-1:]
                    else:
                        cur_seq += gfapy.sequence.rc(gfa.segments[segmap[suc]].sequence)[klen-1:]

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
