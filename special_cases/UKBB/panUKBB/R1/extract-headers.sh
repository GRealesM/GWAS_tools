#!/bin/bash

shopt -s nullglob
array=(*.gz)

for f in ${array[@]};
do
echo -e "$f\t$(zcat $f | head -n1)" >> headerfile.txt
done
