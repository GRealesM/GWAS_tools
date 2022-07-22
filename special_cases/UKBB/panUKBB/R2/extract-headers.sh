#!/bin/bash

shopt -s nullglob
array=(ds/*.gz)

for f in ${array[@]};
do
g=$(echo $f | sed 's/ds\///') # Rm directory from name
echo -e "$g\t$(zcat $f | head -n1)" >> headerfile.txt
done

grep "beta_meta_hq" headerfile.txt | awk '{print $1, "meta_hq"}' >> files_by_header
grep -v "beta_meta_hq" headerfile.txt | grep "beta_meta" | awk '{print $1, "meta_nohq"}' >> files_by_header
grep -vw "beta_meta_hq\|beta_meta" headerfile.txt | grep "beta_EUR" | awk '{print $1, "EUR"}' >> files_by_header
grep -v "beta_meta_hq\|beta_meta\|beta_EUR" headerfile.txt | awk '{print $1, "AFR"}' >> files_by_header
