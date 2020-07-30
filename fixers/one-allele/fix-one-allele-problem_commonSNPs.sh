#!/bin/bash

# Sometimes files come with only one allele (usually, the effect allele). This script is intended to use 1000 genomes data to 
# extract both alleles for all SNPs in the target dataset, and then assign the other allele to REF (in case we know the one allele is the ALT).
# In some cases we have many files with many SNPs in common, since doing the extracting process for the same alleles many times is a waste of time and resources, here's an alternative.

array=(~/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3*legend.gz)
# We define here the initial files to look up for SNPs
initialtargets=(TNFa_Traglia_30134952_1.tsv.gz TNFa_Traglia_30134952_2.tsv.gz TNFa_Traglia_30134952_3.tsv.gz TNFa_Traglia_30134952_4.tsv.gz)


for target in ${initialtargets[@]};
do

echo "Working on $target"
echo Extracting list of SNPs


CHRCOL=$(zcat $target | awk '
 {
   for(i=1;i<=NF;i++) {
     if($i == "CHR")
       printf(i)
   }
   exit 0
 }
 ')

BPCOL=$(zcat $target | awk '
 {
   for(i=1;i<=NF;i++) {
     if($i == "BP")
       printf(i)
   }
   exit 0
 }
 ')

# We first create a CHR/BP list in target to extract those positions from legend files.
# NOTE: Target files should be in hg19 for this approach to work.
zcat $target | awk -v chrcol="$CHRCOL" -v bpcol="$BPCOL" 'NR>1{print $chrcol":"$bpcol}' | sed -e '/^$/d' -e 's/$/:/' >> coords_to_extract.txt

done
# Then, we removed all duplicates from this list
awk '!seen[$0]++' coords_to_extract.txt > tmp.txt && mv tmp.txt coords_to_extract.txt

echo Extracting both alleles from 1000 genomes data...
	for f in ${array[@]};
	do
	  chrnum=$(echo $f | sed -E 's/.*chr([0-9]+).legend.gz/\1/')
	  grep "^$chrnum:" coords_to_extract.txt > tmpchr.txt
	  zcat $f | cut -d" " -f1 | awk -v chrnum="$chrnum" 'BEGIN{OFS=FS=":"}{print chrnum,$2,$1,$3,$4}' | grep -F -f tmpchr.txt  >> alleles.txt
	done


echo Computing the other allele...
Rscript --vanilla fix-alleles_alt.R

rm alleles.txt coords_to_extract.txt tmpchr.txt

