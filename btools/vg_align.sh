#!/bin/bash
if [ "$#" -lt 3 ]; then
    echo "script.sh <graph.gfa> <queries.fastq> <output_folder> [k = 15]"
    exit
fi

graph=$1
queries=$2
ograph=$3/$(basename $graph .gfa)
split=$ograph.split
oqueries=$3/$(basename $queries .fastq)

k=15
if [ "$#" -ge 4 ]; then
    k=$4
fi

vg view -vF $graph > $ograph.vg
vg mod -X 1024 $ograph.vg > $split.vg
vg view -g $split.vg > $split.gfa
vg view -j $split.vg > $split.json
vg index $split.vg -x $split.vg.xg -g $split.vg.gcsa -k $k
vg map -d $split.vg -x $split.vg.xg -g $split.vg.gcsa -f $queries > $oqueries.gam
vg view -aj $oqueries.gam > $oqueries.json
