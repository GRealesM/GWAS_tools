#!/bin/bash

urlarray=($(cat urls_ukbb.txt |tr "\n" " "))
filenamearray=($(cat names_ukbb.txt |tr "\n" " "))

for i in ${!urlarray[@]};
do
echo "Downloading "${filenamearray[$i]}""
wget -O "${filenamearray[$i]}"_Neale_UKBB_1.tsv.gz "${urlarray[$i]}" 
zcat "${filenamearray[$i]}"_Neale_UKBB_1.tsv.gz | head -n1 > header.txt
zcat "${filenamearray[$i]}"_Neale_UKBB_1.tsv.gz | grep -E -f mhc_manifest.txt > body.txt
cat header.txt body.txt | gzip > tmp.tsv.gz && mv tmp.tsv.gz "${filenamearray[$i]}"_Neale_UKBB_1.tsv.gz 
rm header.txt body.txt
echo "Done. File is stored in "${filenamearray[$i]}"_Neale_UKBB_1.tsv.gz"
done

Rscript Preprocess_UKBB.R

