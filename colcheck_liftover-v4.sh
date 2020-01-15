#!/bin/bash

# This code is aimed to check the builds in the files

# Grep using a column in manifest (SNPID) - find SNPs present in the file AND in the manifest
# If there's at least one match, pick first line and check value in BP column
# Compare this value with the manifest file, and check which column matches
# Repeat this process for 10 SNPs (or so), if they agree, set build for overlift.
# This process will be merged with the liftover script, so we identify the build and liftover in one step.


shopt -s nullglob
files=(*.tsv.gz)
#files=(ALR_Ishigaki_doi101101795948_1-hc.tsv.gz)

awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1}' Manifest_build_translator.tsv > rs_manifest.txt
awk 'BEGIN{FS="\t";OFS=":"} NR>1 {print $2,$3}' Manifest_build_translator.tsv > hg18_manifest.txt
awk 'BEGIN{FS="\t";OFS=":"} NR>1 {print $4,$5}' Manifest_build_translator.tsv > hg19_manifest.txt
awk 'BEGIN{FS="\t";OFS=":"} NR>1 {print $6,$7}' Manifest_build_translator.tsv > hg38_manifest.txt


for f in ${files[@]};
do

# Minimum column checkpoint and calculation of beta if necessary
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
	echo $f seem to lack CHR/BP coordinates. | tee -a issues.log
elif [[ "$minSNPID" == 0 ]]; then
	echo $f seem to lack SNPIDs... | tee -a issues.log
elif [[ "$minREF" == 0 || "$minALT" == 0 ]]; then
	echo $f seem to lack REF/ALT labels... | tee -a issues.log
elif [[ "$minP" == 0 ]]; then
	echo $f seem to lack P and SE labels... | tee -a issues.log
elif [[ "$minOR" == 0 && "$minBETA" == 0 ]]; then 
	echo $f seem to lack OR AND BETA... | tee -a issues.log
elif [[  "$minOR" == 1 && "$minBETA" == 0 ]]; then
	echo $f lacks BETA but has OR. BETA will be calculated | tee -a issues.log 
	ORCOL=`zcat $f | awk -F'\t' '
		{
	  		for(i=1;i<=NF;i++) {
	    		if($i == "OR")
	      		printf(i)
	  		}
	  		exit 0
		}
		'`
	zcat $f | awk -v orcol="$ORCOL" 'BEGIN{FS="\t";OFS="\t"} {print $0,log($orcol)}' | sed '1s/-inf/BETA/' > tmp_file_betachecked.tsv
else
	echo $f has all required columns. Excellent!
	zcat $f > tmp_file_betachecked.tsv
fi

### LIFTOVER CHUNK



# NEXT VARIABLES Not necessary if it was defined already

FILEBASENAME=$(echo $f | cut -d- -f1) # Take as file base name everything before the first dash
CHRCOL=`cat tmp_file_betachecked.tsv | awk -F'\t' '
{
  for(i=1;i<=NF;i++) {
    if($i == "CHR")
      printf(i)
  }
  exit 0
}
'`
BPCOL=`cat tmp_file_betachecked.tsv | awk -F'\t' '
{
  for(i=1;i<=NF;i++) {
    if($i == "BP")
      printf(i)
  }
  exit 0
}
'`
SNPIDCOL=`cat tmp_file_betachecked.tsv | awk -F'\t' '
{
  for(i=1;i<=NF;i++) {
    if($i == "SNPID")
      printf(i)
  }
  exit 0
}
'`
# Prepare input for liftover
# Extract relevant columns for target file, and create a BED file containing them. This will be the file fed to the liftover script 
echo "Preparing input for liftover."

cat tmp_file_betachecked.tsv | awk  -v snpidcol="$SNPIDCOL" -v chrcol="$CHRCOL" -v bpcol="$BPCOL" 'BEGIN{FS="\t";OFS="\t"}{print "chr"$chrcol, $bpcol, (($bpcol + 1)), $snpidcol }' | tail -n+2 > ${FILEBASENAME}.bed

# In addition, we'll need to edit the file to make it "mergeable", since merge command will naturally try to merge using the first column in each file. We'll guarantee that the first column in both files to be SNPID.
cp tmp_file_betachecked.tsv tmp_formerging0.tsv

# We reorder the columns only if SNPID is not the first one
if [[ "$SNPIDCOL" -gt 1 ]]; then
	paste <(cut -f${SNPIDCOL} tmp_formerging0.tsv ) <(cut -f1-$((${SNPIDCOL} - 1)),$((${SNPIDCOL} + 1))- tmp_formerging0.tsv) > tmp_formerging1.tsv 
	rm tmp_formerging0.tsv
else
	mv tmp_formerging0.tsv tmp_formerging1.tsv
fi

# Check build
echo Checking build 

# First we'll check if the file is an harmonized file, in which case we'll skip the harmonizing step
harmonized=( $(grep -c "\<hm_CHR\>" tempcolcheck.txt))
if [[ $harmonized -gt 0 ]]; then
	echo "Seems that $f is an harmonized file already, I will skip liftover in this file"
	gzip < tmp_file_betachecked.tsv  > ../04-Liftovered/${FILENAMEBASE}-hg38.tsv.gz
	continue
fi

# In general, we'll use rsids as a proxy to identify. If those aren't available we'll need to figure out which build is by using positions

n_rs=( $(awk -v snpidcol="$SNPIDCOL" 'BEGIN{FS=OFS="\t"}NR>1{print $snpidcol}' tmp_file_betachecked.tsv | grep -c rs))

if [[ "$n_rs" == 0 ]]; then
	echo "Couldn't identify valid rsids (eg. rs429358) in $f. Trying to compare positions instead"
	awk -v chrcol="$CHRCOL" -v bpcol="$BPCOL" 'BEGIN{FS="\t";OFS=":"}NR>1{ print $chrcol, $bpcol}' tmp_file_betachecked.tsv > tmp_file_chrbp.tsv
	matches18=( $(grep -wf hg18_manifest.txt -c tmp_file_chrbp.tsv))
	matches19=( $(grep -wf hg19_manifest.txt -c tmp_file_chrbp.tsv))
	matches38=( $(grep -wf hg38_manifest.txt -c tmp_file_chrbp.tsv))
	if [[ $matches38 -gt $matches19 ]] && [[ $matches38 -gt $matches18 ]];then
		CHOSEN_BUILD=7
	elif [[ $matches19 -gt $matches38 ]] && [[ $matches19 -gt $matches18 ]]; then
		CHOSEN_BUILD=5
	else 	CHOSEN_BUILD=3
	fi
	rm tmp_file_chrbp.tsv
else 
	cat tmp_file_betachecked.tsv | grep -wf rs_manifest.txt | awk -v bpcol="$BPCOL" -v snpidcol="$SNPIDCOL" 'BEGIN{FS="\t";OFS="\t"} {print $snpidcol, $bpcol}' > tmp_rsmatchedfile.txt
	rslist=( $( awk 'BEGIN{FS="\t";OFS="\t"} {print $1}' tmp_rsmatchedfile.txt ) )
	bplist=( $( awk 'BEGIN{FS="\t";OFS="\t"} {print $2}' tmp_rsmatchedfile.txt ) )
	BUILDS=()
	
		for index in ${!rslist[*]};
		do
		
		BUILDS+=( $(grep ${rslist[index]} Manifest_build_translator.tsv | awk -v x="${bplist[index]}" 'BEGIN{OFS=FS="\t"}{for (i=1;i<=NF;i++) if($i == x ) print i }') ) 
		
		done
	
	CHOSEN_BUILD=( $(printf '%d\n' "${BUILDS[@]}" | sort -n | uniq -c | sort -k1,1nr | awk 'NR=1 {print $2; exit}'))
	rm tmp_rsmatchedfile.txt
fi


# Note: 3 = hg18, 5 = hg19, and 7 = hg38
if [ $CHOSEN_BUILD -eq 3 ]; then
	echo $f seems to be in hg18.
	echo Liftovering...
	liftOver ${FILEBASENAME}.bed hg18ToHg38.over.chain.gz ${FILEBASENAME}-lo-output.bed ${FILEBASENAME}-unlifted.bed
elif [ $CHOSEN_BUILD -eq 5 ]; then
	echo $f seems to be in hg19.
	echo Liftovering...
	liftOver ${FILEBASENAME}.bed hg19ToHg38.over.chain.gz ${FILEBASENAME}-lo-output.bed ${FILEBASENAME}-unlifted.bed
elif [ $CHOSEN_BUILD -eq 7 ]; then
	echo $f is in hg38 already, skipping liftover step...
	gzip < tmp_file_betachecked.tsv  > ../04-Liftovered/${FILENAMEBASE}-hg38.tsv.gz
	continue # Remove/rethink this break when the whole pipeline is a single loop and there are more steps afterwards.
else
	echo Sorry, something wrong happened at the build selection stage and I could not identify the build for $f.
	echo $f could not be overlifted | tee -a issues.log
	continue 
fi

awk 'BEGIN{FS="\t";OFS="\t"}{print $4,$1,$2}' ${FILEBASENAME}-lo-output.bed | sed 's/chr//' | sed '1i SNPID\tCHR38\tBP38' > ${FILEBASENAME}-lo-output2.bed
join -a2 -e'NA' -t $'\t' --nocheck-order -o auto ${FILEBASENAME}-lo-output2.bed tmp_formerging1.tsv | gzip  > ../04-Liftovered/${FILEBASENAME}-hg38.tsv
echo $f suscessfully lifted over to hg38 build!

rm  tmp_formerging1.tsv *bed tempcolcheck.txt tmp_file_betachecked.tsv

done
rm rs_manifest.txt hg18_manifest.txt hg19_manifest.txt hg38_manifest.txt
