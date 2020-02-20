#!/bin/bash

shopt -s nullglob
array=(*.gz)

for f in ${array[@]};
do
    FILEBASENAME=$(echo $f | cut -d- -f1); # Take as file base name everything before the first dash

###################
## REHEADING STAGE
###################

echo Working on $FILEBASENAME file
echo Reheading...
    zcat $f | awk -F'\t' 'NR=1{\
	gsub(/\<Log10p\>/,"LOG10P");\
	gsub(/\<_-log10_p-value\>/,"-LOG10P");\
	gsub(/\<effect_allele\>|\<EffectAllele\>|\<A1_effect\>|\<RISK_ALLELE\>|\<EA\>|\<Risk_Allele\>|\<EFFECT_ALLELE\>|\<Alt\>/,"ALT");\
	gsub(/\<OtherAllele\>|\<reference_allele\>|\<OTHER_ALLELE\>|\<other_allele\>|\<A2_other\>|\<NEA\>|\<Ref_Allele\>|\<Ref\>/,"REF");\
# Caution! Sometimes "other_allele" means effect allele, check papers prior to run the script, and pre-rename accordingly.
	gsub(/\<Beta\>|\<beta\>|\<Effect\>|\<effect\>|\<EFFECT\>|\<sebeta_SNP_add\>|\<EFFECT_ALT\>/,"BETA");\
	gsub(/\<Pos\>|\<base_pair_location\>|\<BP\>|\<BP\(hg19\)\>|\<Position\>|\<POS\>|\<pos\>|\<Chr_Position\>|\<bp\>|\<position\>|\<Position\(hg19\)\>|\<POSITION\>|\<bp_hg19\>|\<Coordinate\>|\<chrloc\>/,"BP");\
	gsub(/\<Chr\>|\<chromosome\>|\<Chromosome\>|\<chr\>|\<Chr_ID\>|\<hg18chr\>|\<CHROMOSOME\>/,"CHR");\
	gsub(/\<EMP_Beta\>/,"EMP_BETA");\
	gsub(/\<EMP1\>/,"EMP_P");\
	gsub(/\<EMP_se\>/,"EMP_SE");\
	gsub(/\<hm_effect_allele\>/,"hm_ALT");\
	gsub(/\<hm_beta\>/,"hm_BETA");\
	gsub(/\<hm_pos\>/,"hm_BP");\
	gsub(/\<hm_chrom\>/,"hm_CHR");\
	gsub(/\<hm_odds_ratio\>/,"hm_OR");\
	gsub(/\<hm_other_allele\>/,"hm_REF");\
	gsub(/\<hm_rsid\>/,"hm_SNPID");\
	gsub(/\<n\>/,"N");\
	gsub(/\<odds_ratio\>|\<Odds_ratio\>|\<or\>|\<OddsRatio\>|\<OR\(A1\)\>|\<ORX\>/,"OR");\
	gsub(/\<p_value\>|\<P.value\>|\<pvalue\>|\<P-value\>|\<pval\>|\<p.value\>|\<Pval\>|\<PVALUE\>|\<Pvalue\>|\<P_VALUE\>|\<P-val\>|\<p\>|\<All.p.value\>|\<P_value\>|\<p-value\>|\<GC-adjusted_P_\>/,"P");\
	gsub(/\<standard_error\>|\<StdErr\>|\<stderr\>|\<sebeta_SNP_add\>|\<se\>|\<STDERR\>/,"SE");\
	gsub(/\<Rsq\>/,"RSQ");\
	gsub(/\<id\>|\<variant_id\>|\<MarkerName\>|\<SNP\>|\<rsid\>|\<SNP_Name\>|\<snp\>|\<snpid\>|\<rsID\>|\<\\#SNPID\>|\<rs_number\>|\<RSID\>|\<rs\>|\<db_SNP_RS_ID\/Marker\>|\<dbSNP_RS_ID\>/,"SNPID");\
	gsub(/\<MARKER\>|\<íd\>|\<Chr\:Position\>/,"CHR:BP"); print}' > tmp_reheaded_file.tsv
echo Reheading done. Checking columns...


######################
## COLUMN CHECK STAGE
######################

awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1}' Manifest_build_translator.tsv > rs_manifest.txt
awk 'BEGIN{FS="\t";OFS=":"} NR>1 {print $2,$3}' Manifest_build_translator.tsv > hg18_manifest.txt
awk 'BEGIN{FS="\t";OFS=":"} NR>1 {print $4,$5}' Manifest_build_translator.tsv > hg19_manifest.txt
awk 'BEGIN{FS="\t";OFS=":"} NR>1 {print $6,$7}' Manifest_build_translator.tsv > hg38_manifest.txt

head tmp_reheaded_file.tsv -n1 > tempcolcheck.txt

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
	echo $f seem to lack CHR/BP coordinates.
elif [[ "$minSNPID" == 0 ]]; then
	echo $f seem to lack SNPIDs...
elif [[ "$minREF" == 0 || "$minALT" == 0 ]]; then
	echo $f seem to lack REF/ALT labels...
elif [[ "$minP" == 0 ]]; then
	echo $f seem to lack P and SE labels...
elif [[ "$minOR" == 0 && "$minBETA" == 0 ]]; then 
	echo $f seem to lack OR AND BETA...
elif [[  "$minOR" == 1 && "$minBETA" == 0 ]]; then
	echo $f lacks BETA but has OR. BETA will be calculated 
	ORCOL=`awk -F'\t' '
		{
	  		for(i=1;i<=NF;i++) {
	    		if($i == "OR")
	      		printf(i)
	  		}
	  		exit 0
		}
		' tmp_reheaded_file.tsv`
	awk -v orcol="$ORCOL" 'BEGIN{FS="\t";OFS="\t"} {print $0,log($orcol)}' tmp_reheaded_file.tsv | sed '1s/-inf/BETA/' > tmp_betachecked_file.tsv
else
	echo $f has all required columns. Excellent!
	mv tmp_reheaded_file.tsv tmp_betachecked_file.tsv
fi

##################
## LIFTOVER STAGE
#################

CHRCOL=`awk -F'\t' '
{
  for(i=1;i<=NF;i++) {
    if($i == "CHR")
      printf(i)
  }
  exit 0
}
' tmp_betachecked_file.tsv`
BPCOL=`awk -F'\t' '
{
  for(i=1;i<=NF;i++) {
    if($i == "BP")
      printf(i)
  }
  exit 0
}
' tmp_betachecked_file.tsv`
SNPIDCOL=`awk -F'\t' '
{
  for(i=1;i<=NF;i++) {
    if($i == "SNPID")
      printf(i)
  }
  exit 0
}
' tmp_betachecked_file.tsv`

# Prepare input for liftover
# Extract relevant columns for target file, and create a BED file containing them. This will be the file fed to the liftover script 
echo "Preparing input for liftover."

cat tmp_betachecked_file.tsv | awk  -v snpidcol="$SNPIDCOL" -v chrcol="$CHRCOL" -v bpcol="$BPCOL" 'BEGIN{FS=OFS="\t"}{print "chr"$chrcol, $bpcol, (($bpcol + 1)), $snpidcol }' | tail -n+2 > ${FILEBASENAME}.bed

# In addition, we'll need to edit the file to make it "mergeable", since merge command will naturally try to merge using the first column in each file. We'll guarantee that the first column in both files to be SNPID.
cp tmp_betachecked_file.tsv tmp_formerging0.tsv

# We reorder the columns only if SNPID is not the first one
if [[ "$SNPIDCOL" -gt 1 ]]; then
	paste <(cut -f${SNPIDCOL} tmp_formerging0.tsv ) <(cut -f1-$((${SNPIDCOL} - 1)),$((${SNPIDCOL} + 1))- tmp_formerging0.tsv) > tmp_formerging1.tsv 
	rm tmp_formerging0.tsv
else
	mv tmp_formerging0.tsv tmp_formerging1.tsv
fi

# Check build
echo Checking build 

# In general, we'll use rsids as a proxy to identify. If those aren't available we'll need to figure out which build is by using positions

n_rs=( $(awk -v snpidcol="$SNPIDCOL" 'BEGIN{FS=OFS="\t"}NR>1{print $snpidcol}' tmp_betachecked_file.tsv | grep -cP "[\n\t]rs[0-9]+[\n\t]"))

if [[ "$n_rs" -eq 0 ]]; then
	echo "Couldn't identify valid rsids (eg. rs429358) in $f. Trying to compare positions instead."
	awk -v chrcol="$CHRCOL" -v bpcol="$BPCOL" 'BEGIN{FS="\t";OFS=":"}NR>1{ print $chrcol, $bpcol}' tmp_betachecked_file.tsv > tmp_file_chrbp.tsv
	matches18=( $(grep -wf hg18_manifest.txt -c tmp_file_chrbp.tsv))
	matches19=( $(grep -wf hg19_manifest.txt -c tmp_file_chrbp.tsv))
	matches38=( $(grep -wf hg38_manifest.txt -c tmp_file_chrbp.tsv))
	if [[ $matches38 -gt $matches19 ]] && [[ $matches38 -gt $matches18 ]];then
		CHOSEN_BUILD=7
	elif [[ $matches19 -gt $matches38 ]] && [[ $matches19 -gt $matches18 ]]; then
		CHOSEN_BUILD=5
	else 	CHOSEN_BUILD=3
	fi
	echo "$f had $matches18 with hg18 build, $matches19 with hg19, and $matches38 with hg38."
	rm tmp_file_chrbp.tsv
else 
	cat tmp_betachecked_file.tsv | grep -wf rs_manifest.txt | awk -v bpcol="$BPCOL" -v snpidcol="$SNPIDCOL" 'BEGIN{FS="\t";OFS="\t"} {print $snpidcol, $bpcol}' > tmp_rsmatchedfile.txt
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
	gzip < tmp_betachecked_file.tsv  > ${FILEBASENAME}-hg38.tsv.gz
	continue # Remove/rethink this break when the whole pipeline is a single loop and there are more steps afterwards.
else
	echo Sorry, something wrong happened at the build selection stage and I could not identify the build for $f.
	echo $f could not be overlifted
	continue 
fi

awk 'BEGIN{FS="\t";OFS="\t"}{print $4,$1,$2}' ${FILEBASENAME}-lo-output.bed | sed 's/chr//' | sed '1i SNPID\tCHR38\tBP38' > ${FILEBASENAME}-lo-output2.bed
join -a2 -e'NA' -t $'\t' --nocheck-order -o auto ${FILEBASENAME}-lo-output2.bed tmp_formerging1.tsv | gzip  > ${FILEBASENAME}-hg38.tsv.gz
echo $f suscessfully lifted over to hg38 build!

 rm  tmp_formerging1.tsv *bed tempcolcheck.txt tmp_betachecked_file.tsv

done
 rm rs_manifest.txt hg18_manifest.txt hg19_manifest.txt hg38_manifest.txt

