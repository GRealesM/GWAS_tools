#!/bin/bash

## Backliftover from hg38 to hg19
# Version 2.0
# Background: This script is intended to `back-liftOver` files that have hg38 coordinates (i.e. CHR38, BP38) but lack hg19 coordinates, necessary in some contexts.
#
# Version 2.0 Updates
# - I made the file more generalisable to take one file at a time and run in array batch jobs.

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echo "Compute missing allele frequencies using a selected population from 1000 Genomes project"
   echo
   echo "Syntax: Compute_freqs.sh [[-f file][-p population]] | [-h]"
   echo "options:"
   echo "f     File to be processed"
   echo "h     Print this Help."
   echo
}


while getopts f:h option
do
case "${option}" in
	f) FILE=${OPTARG};;
	h) 
		Help 
		exit;;
	\?) # incorrect option
        	 echo "Error: Invalid option"
	 	Help
         	exit 1;;
esac
done

 FILEBASENAME=$(echo "$FILE" | cut -d. -f1 | cut -d- -f1); # Take as file base name everything before the first dot or first dash
 zcat $FILE > "$FILEBASENAME"_workingfile.tsv

 CHRCOL=$(head -n1 "$FILEBASENAME"_workingfile.tsv | awk -F'\t' '
 {
   for(i=1;i<=NF;i++) {
     if($i == "CHR38")
       printf(i)
   }
   exit 0
 }
 ')
 BPCOL=$(head -n1 "$FILEBASENAME"_workingfile.tsv | awk -F'\t' '
 {
   for(i=1;i<=NF;i++) {
     if($i == "BP38")
       printf(i)
   }
   exit 0
 }
 ')
 
 SNPIDCOL=$(head -n1 "$FILEBASENAME"_workingfile.tsv | awk -F'\t' '
 {
   for(i=1;i<=NF;i++) {
     if($i == "SNPID")
       printf(i)
   }
   exit 0
 }
 ')
 
# I added a remove NA sed command below to avoid problems at liftover stage
 awk -v snpidcol="$SNPIDCOL" -v chrcol="$CHRCOL" -v bpcol="$BPCOL" 'BEGIN{FS=OFS="\t"}{print "chr"$chrcol, $bpcol, (($bpcol + 1)), $snpidcol }' "${FILEBASENAME}"_workingfile.tsv | sed '/NA/d' | tail -n+2 > $FILEBASENAME.bed

 ~/rds/rds-cew54-basis/GWAS_tools/01-Pipeline/liftOver "$FILEBASENAME".bed ~/rds/rds-cew54-basis/GWAS_tools/01-Pipeline/hg38ToHg19.over.chain.gz "$FILEBASENAME"-lo-output.bed "$FILEBASENAME"-unlifted.bed
 awk 'BEGIN{FS=OFS="\t"}{sub("chr", "",$1); print $4,$1,$2}' "$FILEBASENAME"-lo-output.bed | sed '1i SNPID\tCHR19\tBP19' > "$FILEBASENAME"-lo-output2.bed

# Prior to merging, we reorder the columns only if SNPID is not the first one
 if [[ "$SNPIDCOL" -gt 1 ]]; then
 	paste <(cut -f"${SNPIDCOL}" "$FILEBASENAME"_workingfile.tsv ) <(cut -f1-$((SNPIDCOL - 1)),$((SNPIDCOL + 1))- "$FILEBASENAME"_workingfile.tsv) > tmp.tsv && mv tmp.tsv "$FILEBASENAME"_workingfile.tsv 
 fi

 awk -v OFS='\t' 'NR==FNR{a1[$1]=$2; a2[$1]=$3;next};{ if ($1 in a1) print $0, a1[$1], a2[$1]; else print $0, "NA","NA"}' "$FILEBASENAME"-lo-output2.bed "$FILEBASENAME"_workingfile.tsv | gzip  > "$FILEBASENAME"-hg38with19.tsv.gz

 #join -a1 -e'NA' -t $'\t' --nocheck-order -o auto "${FILEBASENAME}"-lo-output2.bed workingfile.tsv | awk '!seen[$0]++' | gzip  > "${FILEBASENAME}"-hg38with19.tsv.gz
 #snpsbeforeliftover=$(echo ""$(cat workingfile.tsv | wc -l)" - 1" | bc)
# snpsafterliftover=$(echo ""$(zcat "${FILEBASENAME}"-hg38with19.tsv.gz | wc -l)" -1" | bc)
# snpsdifference=$(echo "$snpsbeforeliftover" - "$snpsafterliftover" | bc)
# echo ""$f" had "$snpsbeforeliftover" SNPs. After liftover it now has "$snpsafterliftover" ("$snpsdifference" less)."
 echo "$FILE" suscessfully back-lifted over to hg19 build!

 rm "$FILEBASENAME"_workingfile.tsv "$FILEBASENAME"-lo-*.bed "$FILEBASENAME".bed "$FILEBASENAME"-unlifted.bed

