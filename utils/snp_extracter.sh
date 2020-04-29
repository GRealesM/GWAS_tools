#!/bin/bash

shopt -s nullglob
files=(*.tsv.gz)

for f in ${files[@]};
do

SNPIDCOL=`zcat $f | awk -F'\t' '
{
  for(i=1;i<=NF;i++) {
    if($i == "SNPID")
      printf(i)
  }
  exit 0
}
'`
# Capture SNPID column and put it in file 
zcat $f | awk  -v snpidcol="$SNPIDCOL" 'BEGIN{FS="\t";OFS="\t"}{print $snpidcol}' | tail -n+2 >> snp_id_container.txt
done

cat snp_id_container.txt | sort | uniq -c | sort -rn > snp_counter.txt
