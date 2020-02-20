#!/bin/bash

# This code is aimed to check the builds in the files

# Grep using a column in manifest (SNPID) - find SNPs present in the file AND in the manifest
# If there's at least one match, pick first line and check value in BP column
# Compare this value with the manifest file, and check which column matches
# Repeat this process for 10 SNPs (or so), if they agree, set build for overlift.
# This process will be merged with the liftover script, so we identify the build and liftover in one step.

shopt -s nullglob
files=(*.tsv.gz)

awk 'BEGIN{FS="\t";OFS="\t"} NR>1 {print $1}' Manifest_build_translator.tsv > rs_manifest.txt


for f in ${files[@]};
do

# Not necessary if it was defined already

FILEBASENAME=$(echo $f | cut -d- -f1) # Take as file base name everything before the first dash
CHRCOL=`zcat $f | awk -F'\t' '
{
  for(i=1;i<=NF;i++) {
    if($i == "CHR")
      printf(i)
  }
  exit 0
}
'`
BPCOL=`zcat $f | awk -F'\t' '
{
  for(i=1;i<=NF;i++) {
    if($i == "BP")
      printf(i)
  }
  exit 0
}
'`
SNPIDCOL=`zcat $f | awk -F'\t' '
{
  for(i=1;i<=NF;i++) {
    if($i == "SNPID")
      printf(i)
  }
  exit 0
}
'`

# Prepare input for overlift
# Extract relevant columns for target file, and create a BED file containing them. This will be the file fed to the liftover script 
echo "Preparing input for overlift."

zcat $f | awk  -v snpidcol="$SNPIDCOL" -v chrcol="$CHRCOL" -v bpcol="$BPCOL" 'BEGIN{FS="\t";OFS="\t"}{print "chr"$chrcol, $bpcol, (($bpcol + 1)), $snpidcol }' | tail -n+2 > ${FILEBASENAME}.bed

# In addition, we'll need to edit the file to make it "mergeable", since merge command will naturally try to merge using the first column in each file. We'll guarantee that the first column in both files to be SNPID.
zcat $f > tmp_formerging0.tsv

# We reorder the columns only if SNPID is not the first one
if [ "$SNPIDCOL" -gt 1 ]; then
	paste <(cut -f${SNPIDCOL} tmp_formerging0.tsv ) <(cut -f1-$((${SNPIDCOL} - 1)),$((${SNPIDCOL} + 1))- tmp_formerging0.tsv) > tmp_formerging1.tsv 
	rm tmp_formerging0.tsv
else
	mv tmp_formerging0.tsv tmp_formerging1.tsv
fi

# Check build
echo "Checking build"
zcat $f | grep -wf rs_manifest.txt | awk -v bpcol="$BPCOL" -v snpidcol="$SNPIDCOL" 'BEGIN{FS="\t";OFS="\t"} {print $snpidcol, $bpcol}' > tmp_rsmatchedfile.txt
rslist=( $( awk 'BEGIN{FS="\t";OFS="\t"} {print $1}' tmp_rsmatchedfile.txt ) )
bplist=( $( awk 'BEGIN{FS="\t";OFS="\t"} {print $2}' tmp_rsmatchedfile.txt ) )
BUILDS=()

	for index in ${!rslist[*]};
	do
	
	BUILDS+=( $(grep ${rslist[index]} Manifest_build_translator.tsv | awk -v x="${bplist[index]}" 'BEGIN{OFS=FS="\t"}{for (i=1;i<=NF;i++) if($i == x ) print i }') ) 
	
	done

CHOSEN_BUILD=( $(printf '%d\n' "${BUILDS[@]}" | sort -n | uniq -c | sort -k1,1nr | awk 'NR=1 {print $2; exit}'))

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
	zcat $f > ${FILENAMEBASE}-hg38.tsv
	break # Remove/rethink this break when the whole pipeline is a single loop
else
	echo Sorry, something wrong happened at the build selection stage and I could not identify the build for $f.
	echo $f could not be overlifted >> issues.log
	break # This will skip the file for the rest of steps, I think I could keep it 
fi

awk 'BEGIN{FS="\t";OFS="\t"}{print $4,$1,$2}' ${FILEBASENAME}-lo-output.bed | sed 's/chr//' | sed '1i SNPID\tCHR38\tBP38' > ${FILEBASENAME}-lo-output2.bed
join -a2 -e'NA' -t $'\t' --nocheck-order -o auto ${FILEBASENAME}-lo-output2.bed tmp_formerging1.tsv > ${FILEBASENAME}-hg38.tsv
echo $f suscessfully lifted over to hg38 build!



done

rm rs_manifest.txt tmp_rsmatchedfile.txt tmp_formerging1.tsv *bed 
