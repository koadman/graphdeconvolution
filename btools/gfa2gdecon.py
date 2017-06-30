#!/usr/bin/env python3
import gfapy
import sys

depthsfile = open(sys.argv[2])
gfa = gfapy.Gfa.from_file(sys.argv[1])
strains = sys.argv[3]
exp_length = sys.argv[4]

klen = int(gfa.header.kk)

# add basic info
print("G<-"+strains)
print("maxdepth<-200") # fixme
print("expected_length<-"+exp_length)
print("V<-"+str(len(gfa.segments)))
samples = 0

# write out the depths as given in Chris' file
# todo: need to work out some standard way of obtaining these from tigops
depthstr = "depths<-c("
sepc = ""
for line in depthsfile:
    dd = line.rstrip().split(",")
    if samples == 0:
        samples = len(dd)-1
    for d in dd[1:]:
        depthstr += sepc + d
        sepc = ","
depthstr += ")"
print("S<-"+str(samples))
print(depthstr)

# write out unitig lengths
lens = "lengths<-c("
lensep = ""
for seg in gfa.segments:
    lens += lensep + str(len(seg.sequence)-klen)
    lensep = ","
print(lens + ")")

# parse and write out unitig adjacency structure
fromfirst = dict()
fromsecond = dict()
for edge in gfa.edges:
    fo = ""
    rfo = "-"
    if edge.from_orient == "-":
        fo = "-"
        rfo = ""
    to = ""
    rto = "-"
    if edge.from_orient == "-":
        to = "-"
        rto = ""

    # assume integer unitig IDs indexed from 0 -- and convert these to 1-based
    fname = fo+str(int(edge.from_segment.name)+1)
    tname = to+str(int(edge.to_segment.name)+1)
    if fname not in fromfirst:
        fromfirst[fname] = tname
    else:
        fromsecond[fname] = tname
    rfname = rfo+str(int(edge.from_segment.name)+1)
    rtname = rto+str(int(edge.to_segment.name)+1)
    if rtname not in fromfirst:
        fromfirst[rtname] = rfname
    else:
        fromsecond[rtname] = rfname


unifrom = str(len(fromfirst)-len(fromsecond))

print("adj1count<-"+unifrom)
print("adj2count<-"+str(len(fromsecond)))

unifromstr = "adj1source<-c("
unitostr = "adj1dest<-c("
unisep = ""
twofromstr = "adj2source<-c("
twotostr1 = "adj2dest1<-c("
twotostr2 = "adj2dest2<-c("
twosep = ""

for seg in fromfirst:
    if seg in fromsecond:
        twofromstr += twosep + seg
        twotostr1 += twosep + fromfirst[seg]
        twotostr2 += twosep + fromsecond[seg]
        twosep = ","
    else:
        unifromstr += unisep + seg
        unitostr += unisep + fromfirst[seg]
        unisep = ","

print(unifromstr+")")
print(unitostr+")")
print(twofromstr+")")
print(twotostr1+")")
print(twotostr2+")")

#print(line.KC)

