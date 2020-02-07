#!/bin/bash

###################################################
## GWAS Summary Statistics Harmonization Pipeline
##
###################################################

## This script is intended to harmonize GWAS summary statistics to a standard format, and includes the following steps:
## (1) - Detection of file separator and replacement by tabs.
## (2) - Detection of column names and substitution by standard ones.
## (3) - Sanity check of minimum columns, and computation of BETA when missing.
## (3) - Detection of genome build and liftover from hg18 or hg19 to hg38.
## (4) - (In construction) Alignment to reference genome and allele flipping.
##

##################################################
## FILE REQUIREMENTS
##################################################

# This pipeline is meant to be used on GWAS summary statistics with the following characteristics
#
# - Compressed (in GZ format).
# - Fields separated by spaces, tabs, commas, or semicolons.
# - Contains information on:
#	Chromosome (CHR)
#	Position (BP)
#	SNP id (SNPID)
#	Reference allele (REF)*
#	Alternative, or effect allele (ALT)*
#	P-value for association (P)
#	Odds-ratio (OR) or BETA (BETA) 
# 	Standard Error of the log OR (SE) 
# Note: Tags in parentheses are suggestions, since the program will try to infer them from a number of possible tag names
# 	(See dictionary below).
# * We've seen wild heterogeneity in Effect/Reference allele nomenclature (eg. A1/A2 meaning either REF/ALT or ALT/REF,
# 	depending on the study) so I recommend to explicitly rename them prior to running the pipeline.
# - Must be in either GRCh36/hg18, GRCh37/hg19, or GRCh38/hg38 builds.
#

shopt -s nullglob

#############################################
## GLOBAL VARIABLES (Generated once per run)
#############################################

array=(*.tsv.gz) # Initial file list, make sure to have all in tsv.gz format, since chain files are .gz, so including only .gz will try to process those, too :/

## IMD Basis SNP ids and coordinates for the 3 available builds, for build detection

awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1}' Manifest_build_translator.tsv > rs_manifest.txt
awk 'BEGIN{FS="\t";OFS=":"} NR>1 {print $2,$3}' Manifest_build_translator.tsv > hg18_manifest.txt
awk 'BEGIN{FS="\t";OFS=":"} NR>1 {print $4,$5}' Manifest_build_translator.tsv > hg19_manifest.txt
awk 'BEGIN{FS="\t";OFS=":"} NR>1 {print $6,$7}' Manifest_build_translator.tsv > hg38_manifest.txt


#######################
## MAIN LOOP
#######################



for f in ${array[@]};
do

 FILEBASENAME=$(echo $f | cut -d. -f1); # Take as file base name everything before the first dot


#########################
## FILE SEPARATOR STAGE
#########################

 echo Working on $FILEBASENAME file
 echo Checking and changing column separator
 
 n_spaces=$(zcat $f | head -n2 | tail -n1 | tr -d -c ' ' | wc -m)
 n_tabs=$(zcat $f | head -n2 | tail -n1 | tr -d -c '\t' | wc -m)
 n_commas=$(zcat $f | head -n2 | tail -n1 | tr -d -c ',' | wc -m)
 n_semicolons=$(zcat $f | head -n2 | tail -n1 | tr -d -c ';' | wc -m)
 
 if [[ $n_tabs -gt $n_spaces ]] && [[ $n_tabs -gt $n_commas ]] && [[ $n_tabs -gt $n_semicolons ]]; then
 	echo "Separator seems to be tabs. No substitution required. Removing only spaces from header (if any)."
 	zcat $f | sed '1s/ /_/g' > tmp_fs_file.tsv
 elif [[ $n_spaces -gt $n_commas ]] && [[ $n_spaces -gt $n_semicolons ]]; then
 	echo Separator seems to be spaces. Replacing...
 	zcat $f | sed -e 's/^ *//' -e 's/ *$//' -e 's/ \{1,\}/\t/g' > tmp_fs_file.tsv
 elif [[ $n_commas -gt $n_semicolons ]]; then
 	echo Separator seems to be commas. Replacing...	
 	zcat $f | sed -e 's/^,*//' -e 's/,*$//' -e 's/,/\t/g' -e '1s/ /_/g' > tmp_fs_file.tsv
 else 
 	echo Separator seems to be semicolons. Replacing...
 	zcat $f | sed -e 's/^;*//' -e 's/;*$//' -e 's/;/\t/g' -e '1s/ /_/g' > tmp_fs_file.tsv
 fi

####################
## REHEADING STAGE
####################

 echo Reheading...
 	head -n1 tmp_fs_file.tsv | awk -F'\t' '{\
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
 	gsub(/\<íd\>|\<id\>|\<variant_id\>|\<MarkerName\>|\<SNP\>|\<rsid\>|\<SNP_Name\>|\<snp\>|\<snpid\>|\<rsID\>|\<\\#SNPID\>|\<rs_number\>|\<RSID\>|\<rs\>|\<db_SNP_RS_ID\/Marker\>|\<dbSNP_RS_ID\>/,"SNPID");\
 	gsub(/\<MARKER\>|\<íd\>|\<Chr\:Position\>/,"CHR:BP"); print}' > tmp_reheaded_file.tsv
 	tail -n+2 tmp_fs_file.tsv >> tmp_reheaded_file.tsv

## Code below is for files with CHR:BP columns. To be developed and tested.

#chrbp_together=$(head tmp_reheaded_file | grep "CHR:BP" -c tmp_reheaded_file)
#
#if [[ $chrbp_together -gt 0 ]]; then
#	chrbp=$(awk -F'\t' '
#		{
#	  		for(i=1;i<=NF;i++) {
#	    		if($i == "CHR:BP")
#	      		printf(i)
#	  		}
#	  		exit 0
#		}
#		' tmp_reheaded_file.tsv)
#	awk -v chrbpcol="$ORCOL" 'BEGIN{FS=OFS="\t"} {print $0,log($orcol)}' tmp_reheaded_file.tsv | sed '1s/-inf/BETA/' > tmp_betachecked_file.tsv

 rm tmp_fs_file.tsv

 echo Reheading done. Checking columns...


######################
## COLUMN CHECK STAGE
######################

 echo Checking minimum columns.
 
 head tmp_reheaded_file.tsv -n1 > tempcolcheck.txt
 
 # Totally essential. Must be in every file
 minCHR=$(grep -c "\<CHR\>\|\<hm_CHR\>" tempcolcheck.txt)
 minBP=$(grep -c "\<BP\>\|\<hm_BP\>" tempcolcheck.txt)
 minSNPID=$(grep -c "\<SNPID\>\|\<hm_SNPID\>" tempcolcheck.txt)
 minREF=$(grep -c "\<REF\>\|\<hm_REF\>" tempcolcheck.txt)
 minALT=$(grep -c "\<ALT\>\|\<hm_ALT\>" tempcolcheck.txt)
 minP=$(grep -c "\<P\>\|\<hm_P\>\|\<EMP_P\>\|\<SE\>\|\<EMP_SE\>\|\<LOG10P\>\|\<-LOG10P\>\|" tempcolcheck.txt)
 minOR=$(grep -c "\<OR\>\|\<hm_OR\>" tempcolcheck.txt)
 minBETA=$(grep -c "\<hm_BETA\>\|\<BETA\>\|\<EMP_BETA\>" tempcolcheck.txt)
 
 if [[ "$minCHR" == 0 || "$minBP" == 0 ]]; then
 	echo $f seem to lack CHR/BP coordinates. These are essential. Please check your file.
 	echo Continuing with next file...
 	continue
 elif [[ "$minSNPID" == 0 ]]; then
 	echo $f seem to lack SNPIDs. These are essential. Please check your file.
 	echo Continuing with next file...
 	continue
 elif [[ "$minREF" == 0 || "$minALT" == 0 ]]; then
 	echo $f seem to lack REF/ALT labels.These are essential. Please check your file.
 	echo Continuing with next file...
 	continue
 elif [[ "$minP" == 0 ]]; then
 	echo $f seem to lack P and SE labels.These are essential. Please check your file.
 	echo Continuing with next file...
 	continue
 elif [[ "$minOR" == 0 && "$minBETA" == 0 ]]; then 
 	echo $f seem to lack OR AND BETA. At least one of them is essential. Please check your file.
 	echo Continuing with next file...
 	continue
 elif [[  "$minOR" == 1 && "$minBETA" == 0 ]]; then
 	echo $f lacks BETA but has OR. BETA will be calculated 
 	ORCOL=$(awk -F'\t' '
 		{
 	  		for(i=1;i<=NF;i++) {
 	    		if($i == "OR")
 	      		printf(i)
 	  		}
 	  		exit 0
 		}
 		' tmp_reheaded_file.tsv)
 	awk -v orcol="$ORCOL" 'BEGIN{FS="\t";OFS="\t"} {print $0,log($orcol)}' tmp_reheaded_file.tsv | sed '1s/-inf/BETA/' > tmp_betachecked_file.tsv
 else
 	echo $f has all required columns. Excellent!
 	cp tmp_reheaded_file.tsv tmp_betachecked_file.tsv
 fi

###################
## LIFTOVER STAGE
###################

 CHRCOL=$(awk -F'\t' '
 {
   for(i=1;i<=NF;i++) {
     if($i == "CHR")
       printf(i)
   }
   exit 0
 }
 ' tmp_betachecked_file.tsv)
 BPCOL=$(awk -F'\t' '
 {
   for(i=1;i<=NF;i++) {
     if($i == "BP")
       printf(i)
   }
   exit 0
 }
 ' tmp_betachecked_file.tsv)
 SNPIDCOL=$(awk -F'\t' '
 {
   for(i=1;i<=NF;i++) {
     if($i == "SNPID")
       printf(i)
   }
   exit 0
 }
 ' tmp_betachecked_file.tsv)

# Prepare input for liftover
# Extract relevant columns for target file, and create a BED file containing them. This will be the file fed to the liftover script 
 echo Preparing input for liftover...

# A couple annoying things that we should take care of are
# (1) In some files the X chromosome will be encoded as 23. Apparently the almighty liftOver executable is not very happy with this notation
# (2) Likewise, in some files, BP are expressed as scientific notation, which also makes liftOver cringe in disgust.
# so we'll need to add a little workaround
 awk -v chrcol="$CHRCOL" -v bpcol="$BPCOL" 'BEGIN{FS=OFS="\t"}{sub(/23/,"X",$chrcol)}NR>1{$bpcol=sprintf("%i", $bpcol)}1' tmp_betachecked_file.tsv > tmp_Xchrandbpchecked_file.tsv 
 mv tmp_Xchrandbpchecked_file.tsv tmp_betachecked_file.tsv
 
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
# echo Checking build...I will try to use SNPIDs to check for build. 
 echo Checking build...Using cupcake manifest SNPs to compare genomic coordinates. 

# In general, we'll use rsids as a proxy to identify. If those aren't available we'll need to figure out which build is by using positions

# n_rs=( $(awk -v snpidcol="$SNPIDCOL" 'BEGIN{FS=OFS="\t"}NR>1{print $snpidcol}' tmp_betachecked_file.tsv | grep -cP "[\n\t]rs[0-9]+[\n\t]"))
 
# if [[ "$n_rs" -eq 0 ]]; then
#	echo "Couldn't identify valid rsids (eg. rs429358) for $f. Trying to compare positions instead."
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
 	echo "$f had $matches18 with hg18, $matches19 with hg19, and $matches38 with hg38 for manifest SNPs."
 	rm tmp_file_chrbp.tsv
# else 
# 	cat tmp_betachecked_file.tsv | grep -wf rs_manifest.txt | awk -v bpcol="$BPCOL" -v snpidcol="$SNPIDCOL" 'BEGIN{FS="\t";OFS="\t"} {print $snpidcol, $bpcol}' > tmp_rsmatchedfile.txt
# 	rslist=( $( awk 'BEGIN{FS="\t";OFS="\t"} {print $1}' tmp_rsmatchedfile.txt ) )
# 	bplist=( $( awk 'BEGIN{FS="\t";OFS="\t"} {print $2}' tmp_rsmatchedfile.txt ) )
# 	BUILDS=()
# 	for index in ${!rslist[*]};
# 	do
# 		
# 		BUILDS+=( $(grep ${rslist[index]} Manifest_build_translator.tsv | awk -v x="${bplist[index]}" 'BEGIN{OFS=FS="\t"}{for (i=1;i<=NF;i++) if($i == x ) print i }') ) 
# 		
# 	done
# 
# 	CHOSEN_BUILD=( $(printf '%d\n' "${BUILDS[@]}" | sort -n | uniq -c | sort -k1,1nr | awk 'NR=1 {print $2; exit}'))
# 	rm tmp_rsmatchedfile.txt
# fi

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
 
 awk 'BEGIN{FS=OFS="\t"}{sub("chr", "",$1); print $4,$1,$2}' ${FILEBASENAME}-lo-output.bed | sed '1i SNPID\tCHR38\tBP38' > ${FILEBASENAME}-lo-output2.bed
 join -a2 -e'NA' -t $'\t' --nocheck-order -o auto ${FILEBASENAME}-lo-output2.bed tmp_formerging1.tsv | gzip  > ${FILEBASENAME}-hg38.tsv.gz
 echo $f suscessfully lifted over to hg38 build!
 
 rm  tmp_formerging1.tsv *.bed tempcolcheck.txt tmp_betachecked_file.tsv tmp_reheaded_file.tsv



done
rm rs_manifest.txt hg18_manifest.txt hg19_manifest.txt hg38_manifest.txt

