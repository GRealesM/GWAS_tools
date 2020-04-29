#!/bin/bash

# This script is intended to process FinnGen files.
# FinnGen files are many, already in hg38, and have all required columns - so no need for calculations.

# Thus, I just need to reheader them.


shopt -s nullglob
array=(*.tsv.gz)


for f in "${array[@]}";
do
 echo "File $f"
 FILEBASENAME=$(echo "$f" | cut -d. -f1 | cut -d- -f1); # Take as file base name everything before the first dot or first dash
 zcat $f | sed -e '1s/CHR/CHR38/' -e'1s/BP/BP38/' | gzip > tmp.tsv.gz && mv tmp.tsv.gz ${FILEBASENAME}-hg38.tsv.gz
 echo "Done!"
 done
