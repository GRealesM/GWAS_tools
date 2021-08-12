#!/bin/bash

# Let's extract the alleles from .legend.gz files

array=(../../reference/1000GP_Phase3/1000GP_Phase3*legend.gz)
targets=(*.tsv.gz)

SNPIDCOL=$(zcat $targets | awk '
 {
   for(i=1;i<=NF;i++) {
     if($i == "SNP")
       printf(i)
   }
   exit 0
 }
 ')
zcat $targets | awk -v snpidcol="$SNPIDCOL" 'NR>1{print $snpidcol}' | sed 's/$/:/' > SNP_to_extract.txt


for f in ${array[@]};
do
  zcat $f | grep -F -f SNP_to_extract.txt | cut -d" " -f1  >> Alleles.txt
done 

