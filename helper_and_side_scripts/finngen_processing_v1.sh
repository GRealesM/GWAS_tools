#!/bin/bash

# This script is intended to process FinnGen files.
# FinnGen files are many, already in hg38, and have all required columns - so no need for calculations.

# Thus, I just need to reheader them properly.


shopt -s nullglob
array=(*.tsv.gz)


for f in "${array[@]}";
do
 echo "File $f"
 FILEBASENAME=$(echo "$f" | cut -d. -f1 | cut -d- -f1); # Take as file base name everything before the first dot or first dash
 zcat $f | sed -e '1s/#chrom/CHR38/' -e'1s/pos/BP38/' -e '1s/rsids/SNPID/' -e '1s/pval/P/' -e '1s/beta/BETA/' -e '1s/sebeta/SE' | gzip > ${FILEBASENAME}-hg38.tsv.gz
 echo "Done!"
 done
