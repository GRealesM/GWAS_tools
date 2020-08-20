#!/bin/bash

shopt -s nullglob
#array=(*.gz)
# New array starting from 4th entry
array=(*hg38.tsv.gz)


for f in ${array[@]};
do 

echo $f
FILEBASENAME=$(echo $f | cut -d- -f1); # Take as file base name everything before the first dash

zcat $f | head -n1 | sed -e 's/\bCHR\b/CHR38/' -e 's/\bBP\b/BP38/' > ${FILEBASENAME}-corr.tsv
zcat $f | tail -n+2 >> ${FILEBASENAME}-corr.tsv
gzip ${FILEBASENAME}-corr.tsv

done
