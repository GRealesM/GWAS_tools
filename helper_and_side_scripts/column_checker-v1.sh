#!/bin/bash

# This script will check for issues in files, regarding column names
shopt -s nullglob
files=(*.tsv.gz)

for f in ${files[@]};
do

zcat $f | head -n1 > tempcolcheck.txt

# Totally essential. Must be in every file
minCHR=`grep -c "\<CHR\>\|\<hm_CHR\>" tempcolcheck.txt`
minBP=`grep -c "\<BP\>\|\<hm_BP\>" tempcolcheck.txt`
minSNPID=`grep -c "\<SNPID\>\|\<hm_SNPID\>" tempcolcheck.txt`
minREF=`grep -c "\<REF\>\|\<hm_REF\>" tempcolcheck.txt`
minALT=`grep -c "\<ALT\>\|\<hm_ALT\>" tempcolcheck.txt`
minP=`grep -c "\<P\>\|\<hm_P\>\|\<EMP_P\>\|\<SE\>\|\<EMP_SE\>\|\<LOG10P\>\|\<-LOG10P\>\|" tempcolcheck.txt`
minOR=`grep -c "\<OR\>\|\<hm_OR\>" tempcolcheck.txt`
minBETA=`grep -c "\<hm_BETA\>\|\<BETA\>\|\<EMP_BETA\>" tempcolcheck.txt`

if [[ "$minCHR" == 0 || "$minBP" == 0 ]]; then
	echo "$f seem to lack CHR/BP coordinates."
elif [[ "$minSNPID" == 0 ]]; then
	echo "$f seem to lack SNPIDs..."
elif [[ "$minREF" == 0 || "$minALT" == 0 ]]; then
	echo "$f seem to lack REF/ALT labels..."
elif [[ "$minP" == 0 ]]; then
	echo "$f seem to lack P and SE labels..."
elif [[ "$minOR" == 0 && "$minBETA" == 0 ]]; then 
	echo "$f seem to lack OR AND BETA..."
elif [[  "$minOR" == 1 && "$minBETA" == 0 ]]; then
	echo "$f lacks BETA but has OR. BETA should be calculated"
	# Insert the following code when appropriate
	#	ORCOL=`zcat $f | awk -F'\t' '
	#	{
	#  		for(i=1;i<=NF;i++) {
	#    		if($i == "OR")
	#      		printf(i)
	#  		}
	#  		exit 0
	#	}
	#	'`
	# zcat $f | awk -v orcol="$ORCOL" 'BEGIN{FS="\t";OFS="\t"} {print $0,log($orcol)}' | sed '1s/-inf/BETA/' > DUMMYFILENAME
fi

done

rm tempcolcheck.txt
echo "Done analyzing files!"

