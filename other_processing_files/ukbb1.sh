#!/bin/bash

## Downloading all datasets from UKBB and process them using only the SNPs in the manifests.

# Define the appropriate variables
urlarray=($(cat urls_ukbb.txt |tr "\n" " "))
filenamearray=($(cat names_ukbb.txt |tr "\n" " "))

for i in ${!urlarray[@]};
do
echo "Downloading "${filenamearray[$i]}""
wget -O "${filenamearray[$i]}"_Neale_UKBB_1.tsv.gz "${urlarray[$i]}" 
zcat "${filenamearray[$i]}"_Neale_UKBB_1.tsv.gz | head -n1 > header.txt
zcat "${filenamearray[$i]}"_Neale_UKBB_1.tsv.gz | grep -f key19.tsv > body.txt
cat header.txt body.txt > tmp.tsv && mv tmp.tsv "${filenamearray[$i]}"_Neale_UKBB_1.tsv 
rm header.txt body.txt "${filenamearray[$i]}"_Neale_UKBB_1.tsv.gz
echo "Done. File is stored in "${filenamearray[$i]}"_Neale_UKBB_1.tsv"
done
