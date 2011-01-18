#!/bin/bash
CIRCOS=/users/wallen/src/circos-0.52/bin/circos

python write_chr_files.py
cd chrs
for i in `ls`
do
    $CIRCOS -conf $i &
done
cd ..
