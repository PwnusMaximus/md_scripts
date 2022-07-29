#!/bin/bash

#This script will take n number of stacking data files
#with the left column being the stacking interaction and the right
#being the stacking percentage
#

hdr() { awk 'FNR==1{ print "\0", FILENAME }1' "$1"; }

join -a1 -a2 -o auto <(hdr *.dat) <(hdr *.dat) >join.tmp

for file in FINAL*; do
     join -a1 -a2 -o auto join.tmp <(hdr "$file") >join.tmp.1
     mv join.tmp.1 join.tmp
  done
tr -d '\0' <join.tmp >final.file

echo "Done!" 
echo "Files sorted and merged on column #1"
echo "import 'final.file' into excel for analysis"