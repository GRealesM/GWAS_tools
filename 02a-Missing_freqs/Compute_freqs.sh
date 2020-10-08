#!/bin/bash

# Compute missing allele frequencies using a selected population from 1000 Genomes project
# Version 1.0
#
# This script is intended to extract the allele frequencies of GWAS summary statistics datatests from 1000 genomes data.
# It requires 2 pieces of data, -f the file to be processed, and -p the 1000 genomes population from which to compute the frequencies {ACB,ASW,BEB,CDX,CEU,CHB,CHS,CLM,ESN,FIN,GBR,GIH,GWD,IBS,ITU,JPT,KHV,LWK,MSL,MXL,PEL,PJL,PUR,STU,TSI,YRI}.
# IMPORTANT NOTE 1: Files should be in hg19, like 1000 genomes Phase III files.
# IMPORTANT NOTE 2: Genomic coordinates (hg19) should be denoted by CHR/BP headers.

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
   echo "p     1000 Genomes population code (ACB,ASW,BEB,CDX,CEU,CHB,CHS,CLM,ESN,FIN,GBR,GIH,GWD,IBS,ITU,JPT,KHV,LWK,MSL,MXL,PEL,PJL,PUR,STU,TSI,YRI). Please choose only one."
   echo "h     Print this Help."
   echo
}


while getopts f:p:h option
do
case "${option}" in
	f) FILE=${OPTARG};;
	p) POP=${OPTARG};;
	h) 
		Help 
		exit;;
	\?) # incorrect option
        	 echo "Error: Invalid option"
	 	Help
         	exit 1;;
esac
done

# We'll first extract the SNPs chromosome positions to use them as a query. In this case, all files have the same SNP number, so we can extract this from the first one.
# To make it more column-independent, we'll create variables to detect which columns do we need.

CHRCOL=$(zcat $FILE | awk -F'\t' ' {for(i=1;i<=NF;i++) { if($i == "CHR19") printf(i) } exit 0 }')
BPCOL=$(zcat $FILE | awk -F'\t' '{for(i=1;i<=NF;i++) { if($i == "BP19") printf(i) } exit 0}')
FILENAME=$(echo $FILE | sed -e 's/.tsv.gz//' -e 's/-.*//') 

zcat $FILE | awk -F"\t" -v chrcol="$CHRCOL" -v bpcol="$BPCOL" 'NR>1{print $chrcol,$bpcol,$bpcol,$bpcol}' > hg19snps.txt

# Now we'll get the sample IDs for our population of interest. For that we'll use the integrated_call_samples_v3.20130502.ALL.panel file.
# Please note that this panel makes reference to phase 3 1000 genomes (in hg19). I just couldn't copy it to the appropriate reference folder.
if test -f "~/rds/rds-cew54-basis/GWAS_tools/02a-Missing_freqs/"$POP"_samples.txt"; then
    echo ""$POP"_samples.txt exists."
else grep "$POP" ~/rds/rds-cew54-basis/95-1000genomes/reference_hg19/integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1, $1}' > ~/rds/rds-cew54-basis/GWAS_tools/02a-Missing_freqs/"$POP"_samples.txt
fi 
# We create a header file to host our frequencies.
echo -e "CHR19\tBP19\tSNPID\tREF\tALT\tALT_FREQ\tOBS_CT" > Freqs_"$FILENAME".txt

# We extract frequencies and append to our Freqs file
for chr in {1..22};
do
	echo "Exctracting from chr$chr..."
	grep -P "^"$chr" " hg19snps.txt > tmpsnps.txt
	plink2 --bfile ~/rds/rds-cew54-basis/95-1000genomes/reference_hg19/chr$chr --keep ~/rds/rds-cew54-basis/GWAS_tools/02a-Missing_freqs/"$POP"_samples.txt --extract range tmpsnps.txt --freq cols=+pos --out temp
	tail -n+2 temp.afreq >> Freqs_"$FILENAME".txt
done


Rscript --vanilla ~/rds/rds-cew54-basis/GWAS_tools/02a-Missing_freqs/Mergeback.R $FILE Freqs_"$FILENAME".txt
rm temp* hg19snps.txt tmpsnps.txt Freqs_*

echo "Done!"



