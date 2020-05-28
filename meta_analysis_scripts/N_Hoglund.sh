#!/bin/bash

hogfiles=(../02-Processed/*Hoglund*)

for f in ${hogfiles[@]}; 
do
	echo "Working on $f"
	NCOL=$(zcat $f | head -n1 | awk -F'\t' '
 		{
 	  		for(i=1;i<=NF;i++) {
 	    		if($i == "N")
 	      		printf(i)
 	  		}
 	  		exit 0
 		}
 		')
	N=($(zcat $f | cut -d$'\t' -f"$NCOL" | head -n2 | tail -n+2))
	echo "$f $N" >> N_Hoglund.txt
done
