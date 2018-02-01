#!/usr/bin/env python3
import json
import sys

if len(sys.argv) < 3:
    print("Usage: %s <split_graph.json> <mapping.json>" % sys.argv[0], file=sys.stderr)
    exit(1)

# set up a reverse complement table
fwd_chars = "ACGTacgt"
rev_chars = "TGCAtgca"
rctab = str.maketrans(fwd_chars,rev_chars)

# TODO: this approach of slurping in the whole graph might be too memory hungry
print("Reading assembly graph JSON from "+sys.argv[1], file=sys.stderr)
asm_file = open(sys.argv[1])
asm_json = json.load(asm_file)
seq_id_map = {}
for seq in asm_json["node"]:
    seq_id_map[seq["id"]]=seq["sequence"]

print("Reading mapped reads from " + sys.argv[2], file=sys.stderr)
json_file = open(sys.argv[2])

seq_id_bins = {}
bin_id=-1
for line in json_file:
    bin_id += 1
    gam = json.loads(line)
    name = gam["name"]
    print("Processing sequence " + name, file=sys.stderr)
    if "path" in gam:
        for mapping in gam["path"]["mapping"]:
            if "position" in mapping:
                node_id = mapping["position"]["node_id"]
                if not node_id in seq_id_bins: seq_id_bins[node_id]=set()
                seq_id_bins[node_id].add(bin_id)
            print("Sequence appears to be unaligned", file=sys.stderr)
    else:
        print("Entry for %s sequence had no 'path' attribute" % name, file=sys.stderr)

for node in seq_id_bins:
    node_bins = str(node)
    for b in seq_id_bins[node]:
        node_bins += "\t"+str(b)
    print(node_bins)
