#!/bin/bash

# Let's extract the alleles from .legend.gz files

array=(../../reference/1000GP_Phase3/1000GP_Phase3*legend.gz)

for f in ${array[@]};
do
  zcat $f | grep -F -f SNP_LopezIsac.txt | cut -d" " -f1  >> tomerge_1kg_LopezIsac_alleles.txt
done 

