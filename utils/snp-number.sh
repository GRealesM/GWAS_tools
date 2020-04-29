#!/bin/bash

shopt -s nullglob
array=(*.gz)

for f in ${array[@]};
do
numberoflines=`zcat $f |wc -l`
res=`echo $numberoflines - 1 | bc` 
echo -e "$f\t$res" >> snp-number.txt
done
