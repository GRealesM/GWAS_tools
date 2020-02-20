#!/bin/bash

files=(CRD_Liu_26192919_3-hc.tsv.gz T2D_DIAGRAM_24509480_1-hc.tsv.gz)
builds=(37 36)

# Checking if the list of builds is the same length as the list of files, this probably should be checked at the beginning!
if [ ${#files[@]} != ${#builds[@]} ]; then
	echo "The number of files is ${#files[@]} but you provided ${#builds[@]} builds. Please check! Exiting..."
fi
 

for index in ${!files[*]};
do
# Create filename variable if not present already
FILEBASENAME=$(echo ${files[index]} | cut -d- -f1) # Take as file base name everything before the first dash
# Create variables for CHR, BP, and SNPID column names, as they'll be used later on
CHRCOL=`zcat ${files[index]} | awk -F'\t' '
{
  for(i=1;i<=NF;i++) {
    if($i == "CHR")
      printf(i)
  }
  exit 0
}
'`

BPCOL=`zcat ${files[index]} | awk -F'\t' '
{
  for(i=1;i<=NF;i++) {
    if($i == "BP")
      printf(i)
  }
  exit 0
}
'`

SNPIDCOL=`zcat ${files[index]} | awk -F'\t' '
{
  for(i=1;i<=NF;i++) {
    if($i == "SNPID")
      printf(i)
  }
  exit 0
}
'`
# Prepare input for overlift
zcat ${files[index]} | awk  -v snpidcol="$SNPIDCOL" -v chrcol="$CHRCOL" -v bpcol="$BPCOL" 'BEGIN{FS="\t";OFS="\t"}{print "chr"$chrcol, $bpcol, (($bpcol + 1)), $snpidcol }' | tail -n+2 > ${FILEBASENAME}.bed
zcat ${files[index]} > temp.tsv

if [ "$SNPIDCOL" -gt 1 ]; then
	paste <(cut -f${SNPIDCOL} temp.tsv ) <(cut -f1-$((${SNPIDCOL} - 1)),$((${SNPIDCOL} + 1))- temp.tsv) > temp_for_merging.tsv 
else
	cp temp.tsv temp_for_merging.tsv
fi

# Check builds
if [ ${builds[index]} -eq 36 ]; then
	liftOver ${FILEBASENAME}.bed hg18ToHg38.over.chain.gz ${FILEBASENAME}-lo-output.bed ${FILEBASENAME}-unlifted.bed
elif [ ${builds[index]} -eq 37 ]; then
	liftOver ${FILEBASENAME}.bed hg19ToHg38.over.chain.gz ${FILEBASENAME}-lo-output.bed ${FILEBASENAME}-unlifted.bed
elif [ ${builds[index]} -eq 38 ]; then
	echo File ${files[index]} is in hg38 already, skipping liftover step...
	zcat ${files[index]} > ${FILENAMEBASE}-hg38.tsv
	break # Remove/rethink this break when the whole pipeline is a single loop
else
	echo Sorry, I could not identify ${build[index]} as a valid build. Please check.
	echo ${files[index]} could not be overlifted >> issues.log
	break # This will skip the file for the rest of steps, I think I could keep it 
fi

awk 'BEGIN{FS="\t";OFS="\t"}{print $4,$1,$2}' ${FILEBASENAME}-lo-output.bed | sed 's/chr//' | sed '1i SNPID\tCHR38\tBP38' > ${FILEBASENAME}-lo-output2.bed
join -a2 -e'NA' -t $'\t' --nocheck-order -o auto ${FILEBASENAME}-lo-output2.bed temp_for_merging.tsv > ${FILEBASENAME}-hg38.tsv
echo ${FILEBASENAME} suscessfully lifted over to hg38 build!

done

