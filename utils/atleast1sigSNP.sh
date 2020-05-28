#!/bin/bash

# This script will select all non-metaanalysed files and check if they have at least one significant SNP (P<5e-8)

#shopt -s nullglob

#nonmeta=($(ls | grep "Ahola\|Enroth\|Hoglund"))
nonmeta=(*hg38*)


for f in ${nonmeta[@]};
do
	echo "Working on $f"
	PCOL=$(zcat $f | head -n1 | awk -F'\t' '
 		{
 	  		for(i=1;i<=NF;i++) {
 	    		if($i == "P")
 	      		printf(i)
 	  		}
 	  		exit 0
 		}
 		')
	mostsnp=($(zcat $f | cut -d$'\t' -f"$PCOL" | tail -n+2 | sort -g | head -n1))	
	threshold=5.0E-8
	if (( $(echo "$mostsnp $threshold" | awk '{print ($1 > $2)}') )); then
		echo "$f has no significant SNPs"
		echo $f >> nonsigsnps_datasets.txt;
	else
		echo "$f has at least one significant SNP"
		echo $f >> sigsnps_datasets.txt;
	fi
done
