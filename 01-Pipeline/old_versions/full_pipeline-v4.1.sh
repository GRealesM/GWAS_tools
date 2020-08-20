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



for f in "${array[@]}";
do

 FILEBASENAME=$(echo "$f" | cut -d. -f1); # Take as file base name everything before the first dot


#########################
## FILE SEPARATOR STAGE
#########################

 echo Working on "$FILEBASENAME" file
 echo Checking and changing column separator
 
 n_spaces=$(zcat "$f" | head -n2 | tail -n1 | tr -d -c ' ' | wc -m)
 n_tabs=$(zcat "$f" | head -n2 | tail -n1 | tr -d -c '\t' | wc -m)
 n_commas=$(zcat "$f" | head -n2 | tail -n1 | tr -d -c ',' | wc -m)
 n_semicolons=$(zcat "$f" | head -n2 | tail -n1 | tr -d -c ';' | wc -m)
 
 if [[ $n_tabs -gt $n_spaces ]] && [[ $n_tabs -gt $n_commas ]] && [[ $n_tabs -gt $n_semicolons ]]; then
 	echo "Separator seems to be tabs. No substitution required. Removing only spaces from header (if any)."
 	zcat "$f" | sed '1s/ /_/g' > tmp_fs_file.tsv
 elif [[ $n_spaces -gt $n_commas ]] && [[ $n_spaces -gt $n_semicolons ]]; then
 	echo Separator seems to be spaces. Replacing...
 	zcat "$f" | sed -e 's/^ *//' -e 's/ *$//' -e 's/ \{1,\}/\t/g' > tmp_fs_file.tsv
 elif [[ $n_commas -gt $n_semicolons ]]; then
 	echo Separator seems to be commas. Replacing...	
 	zcat "$f" | sed -e 's/^,*//' -e 's/,*$//' -e 's/,/\t/g' -e '1s/ /_/g' > tmp_fs_file.tsv
 else 
 	echo Separator seems to be semicolons. Replacing...
 	zcat "$f" | sed -e 's/^;*//' -e 's/;*$//' -e 's/;/\t/g' -e '1s/ /_/g' > tmp_fs_file.tsv
 fi

####################
## REHEADING STAGE
####################


### NOTE: I TEMPORARILY ADDED A1 to be replaced by ALT and A2 by REF for Enroth files processing. This may not be the case for other files, so remember to remove them later!


 echo Reheading...
 	head -n1 tmp_fs_file.tsv | awk -F'\t' '{
 	gsub(/\<Log10p\>/,"LOG10P");
 	gsub(/\<_-log10_p-value\>/,"-LOG10P");
 	gsub(/\<effect_allele\>|\<Effect_Allele\>|\<EffectAllele\>|\<A1_effect\>|\<RISK_ALLELE\>|\<EA\>|\<Risk_Allele\>|\<EFFECT_ALLELE\>|\<Alt\>|\<A1\>/,"ALT");
 	gsub(/\<OtherAllele\>|\<reference_allele\>|\<Ref_Allele\>|\<OTHER_ALLELE\>|\<other_allele\>|\<A2_other\>|\<NEA\>|\<Ref_Allele\>|\<Ref\>|\<A2\>/,"REF");
 # Caution! Sometimes "other_allele" means effect allele, check papers prior to run the script, and pre-rename accordingly.
 	gsub(/\<Beta\>|\<beta\>|\<Effect\>|\<effect\>|\<EFFECT\>|\<sebeta_SNP_add\>|\<EFFECT_ALT\>/,"BETA");
 	gsub(/\<Pos\>|\<base_pair_location\>|\<BP\>|\<BP\(hg19\)\>|\<Position\>|\<POS\>|\<pos\>|\<Chr_Position\>|\<bp\>|\<position\>|\<Position\(hg19\)\>|\<POSITION\>|\<bp_hg19\>|\<Coordinate\>|\<chrloc\>/,"BP");
 	gsub(/\<Chr\>|\<chromosome\>|\<Chromosome\>|\<chr\>|\<Chr_ID\>|\<hg18chr\>|\<CHROMOSOME\>|\<chrom\>/,"CHR");
 	gsub(/\<EMP_Beta\>/,"EMP_BETA");
 	gsub(/\<EMP1\>/,"EMP_P");
 	gsub(/\<EMP_se\>/,"EMP_SE");
 	gsub(/\<hm_effect_allele\>/,"hm_ALT");
 	gsub(/\<hm_beta\>/,"hm_BETA");
 	gsub(/\<hm_pos\>/,"hm_BP");
 	gsub(/\<hm_chrom\>/,"hm_CHR");
 	gsub(/\<hm_odds_ratio\>/,"hm_OR");
 	gsub(/\<hm_other_allele\>/,"hm_REF");
 	gsub(/\<hm_rsid\>/,"hm_SNPID");
 	gsub(/\<n\>/,"N");
 	gsub(/\<odds_ratio\>|\<Odds_ratio\>|\<or\>|\<OddsRatio\>|\<OR\(A1\)\>|\<ORX\>/,"OR");
 	gsub(/\<p_value\>|\<P.value\>|\<pvalue\>|\<P-value\>|\<pval\>|\<p.value\>|\<Pval\>|\<PVALUE\>|\<Pvalue\>|\<P_VALUE\>|\<P-val\>|\<p\>|\<All.p.value\>|\<P_value\>|\<p-value\>|\<GC-adjusted_P_\>/,"P");
 	gsub(/\<standard_error\>|\<StdErr\>|\<stderr\>|\<sebeta_SNP_add\>|\<se\>|\<STDERR\>/,"SE");
 	gsub(/\<Rsq\>/,"RSQ");
 	gsub(/\<íd\>|\<id\>|\<variant_id\>|\<MarkerName\>|\<SNP\>|\<rsid\>|\<SNP_Name\>|\<snp\>|\<snpid\>|\<SNP_ID\>|\<rsID\>|\<\\#SNPID\>|\<rs_number\>|\<RSID\>|\<rs\>|\<db_SNP_RS_ID\/Marker\>|\<dbSNP_RS_ID\>/,"SNPID");
 	gsub(/\<MARKER\>|\<íd\>|\<Chr\:Position\>/,"CHR:BP");
	gsub(/\<Zscore\>|\<ZSCORE\>/,"Z");print}' > tmp_reheaded_file.tsv
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

 echo Reheading done.


######################
## COLUMN CHECK STAGE
######################

 echo Checking minimum columns.
 
 head tmp_reheaded_file.tsv -n1 > tempcolcheck.txt
 
 # Totally essential. Must be in every file
 minCHR=$(grep -c "\<CHR\>" tempcolcheck.txt)
 minBP=$(grep -c "\<BP\>" tempcolcheck.txt)
 minSNPID=$(grep -c "\<SNPID\>" tempcolcheck.txt)
 minREF=$(grep -c "\<REF\>" tempcolcheck.txt)
 minALT=$(grep -c "\<ALT\>" tempcolcheck.txt)
 minP=$(grep -c "\<P\>" tempcolcheck.txt)
 minOR=$(grep -c "\<OR\>" tempcolcheck.txt)
 minBETA=$(grep -c "\<BETA\>" tempcolcheck.txt)
 minSE=$(grep -c "\<SE\>" tempcolcheck.txt)
 minZ=$(grep -c "\<Z\>" tempcolcheck.txt)

 
 # In addition, we'll check that essential estimates are not NA (Note: This will only work if the code for missing data is "NA")
 # We define a threshold of acceptance of how many NAs we can accept in each column. We start by 50%, but we can of course reduce this.
 maxNA=0.5

if [[ "$minCHR" == 0 || "$minBP" == 0 ]]; then
 	echo "$f" seem to lack CHR/BP coordinates. These are essential. Please check your file.
 	echo Continuing with next file...
 	continue
 elif [[ "$minSNPID" == 0 ]]; then
 	echo "$f" seem to lack SNPIDs. These are essential. Please check your file.
 	echo Continuing with next file...
 	continue
	
 elif [[ "$minREF" == 0 || "$minALT" == 0 ]]; then
 	echo "$f" seem to lack REF/ALT labels. These are essential. Please check your file.
 	echo Continuing with next file...
 	continue
	
 if [[ "$minP" == 0 ]]; then
 	if [[ "$minZ" == 1 ]]; then
		echo "$f" seems to lack P column, and although it has Z, the functionality for calculating P from Z scores has not been implemented yet, sorry!
		echo Continuing with next file...
		continue
	else 
 		echo "$f" seems to lack P column. This column is essential. Please check your file.
 		echo Continuing with next file...
	 	continue
	fi
  fi
 else # This works for "P" labels only. I'll extend this for more forms of P in the future
  	PCOL=$(awk -F'\t' '
 		{
 	  		for(i=1;i<=NF;i++) {
 	    		if($i == "P")
 	      		printf(i)
 	  		}
 	  		exit 0
 		}
 		' tmp_reheaded_file.tsv)
	NAP=$(echo "scale=2;$(awk -v pcol="$PCOL" 'BEGIN{FS=OFS="\t"}NR>1{print $pcol}' tmp_reheaded_file.tsv | grep -c NA) / ($(cat tmp_reheaded_file.tsv | wc -l)-1)" | bc)
	if [[ 1 -eq $(echo  "$NAP > $maxNA" | bc -l) ]]; then
		echo P column seems to have more NAs than what is acceptable. Please check your file.
		echo Continuing with next file...
		continue
	fi
 fi

 if [[ "$minOR" == 0 && "$minBETA" == 0 ]]; then 
 	echo "$f" seem to lack OR AND BETA. At least one of them is essential. Please check your file.
 	echo Continuing with next file...
 	continue
	fi
 if [[ "$minBETA" == 1 ]]; then
 	BETACOL=$(awk -F'\t' '
 		{
 	  		for(i=1;i<=NF;i++) {
 	    		if($i == "BETA")
 	      		printf(i)
 	  		}
 	  		exit 0
 		}
 		' tmp_reheaded_file.tsv)
	NABETA=$(echo "scale=2;$(awk -v betacol="$BETACOL" 'BEGIN{FS=OFS="\t"}NR>1{print $betacol}' tmp_reheaded_file.tsv | grep -c NA) / ($(cat tmp_reheaded_file.tsv | wc -l)-1)" | bc)
	if [[ 1 -eq $(echo  "$NABETA > $maxNA" | bc -l) ]]; then
		echo BETA column seems to have more NAs than what is acceptable. Please check your file.
		echo Continuing with next file...
		continue
	fi	
fi

# End of first stage. Up to here we have only checked if the veryminimum columns required are present. Some of the columns (ie BETA and SE) can be calculated from available data
mv tmp_reheaded_file.tsv tmp_schecked.tsv

# We can calculate BETA from OR (if OR column exists and is not all NA).
if [[  "$minOR" == 1 && "$minBETA" == 0 ]]; then	
 	ORCOL=$(awk -F'\t' '
 		{
 	  		for(i=1;i<=NF;i++) {
 	    		if($i == "OR")
 	      		printf(i)
 	  		}
 	  		exit 0
 		}
 		' tmp_schecked.tsv)
	NAOR=$(echo "scale=2;$(awk -v orcol="$ORCOL" 'BEGIN{FS=OFS="\t"}NR>1{print $orcol}' tmp_schecked.tsv | grep -c NA) / ($(cat tmp_schecked.tsv | wc -l)-1)" | bc)
	if  [ 1 -eq $(echo  "$NAOR > $maxNA" | bc -l) ]; then
		echo OR column seems to have more NAs than what is acceptable. Please check your file.
		echo Continuing with next file...
		continue	
	else 	
		echo "$f" lacks BETA but has OR. BETA will be calculated 
 		awk -v orcol="$ORCOL" 'BEGIN{FS=OFS="\t"} {print $0,log($orcol)}' tmp_schecked.tsv | sed -e 's///' -e '1s/-inf/BETA/' > tmp.tsv && mv tmp.tsv tmp_schecked.tsv
	fi
fi

# We can calculate SE if we have BETA and Z or P values
if [[ "$minSE" == 0 ]]; then
	echo "$f" seems to lack SE. It will be calculated using BETA and P values. This will calculate Z scores as a side effect too.
	Rscript SE_Calc_pipeline.R
	mv tmp.tsv tmp_schecked.tsv
fi

# A couple annoying things that we should take care of are
# (1) In some files the X chromosome will be encoded as 23. Apparently the almighty liftOver executable is not very happy with this notation
# (2) Likewise, in some files, BP are expressed as scientific notation, which also makes liftOver cringe in disgust.
# (3) In addition, in some files CHR is presented as chrX, which will likely upset liftOver too. We need to remove it.
# so we'll need to add a little workaround
 CHRCOL=$(awk -F'\t' '
 {
   for(i=1;i<=NF;i++) {
     if($i == "CHR")
       printf(i)
   }
   exit 0
 }
 ' tmp_schecked.tsv)
 BPCOL=$(awk -F'\t' '
 {
   for(i=1;i<=NF;i++) {
     if($i == "BP")
       printf(i)
   }
   exit 0
 }
 ' tmp_schecked.tsv)

awk -v chrcol="$CHRCOL" -v bpcol="$BPCOL" 'BEGIN{FS=OFS="\t"}{sub(/23/,"X",$chrcol); sub("chr","",$chrcol)}NR>1{$bpcol=sprintf("%i", $bpcol)}1' tmp_schecked.tsv > tmp.tsv && mv tmp.tsv tmp_schecked.tsv

echo Column sanity check OK. 
	

###################
## LIFTOVER STAGE
###################

 SNPIDCOL=$(awk -F'\t' '
 {
   for(i=1;i<=NF;i++) {
     if($i == "SNPID")
       printf(i)
   }
   exit 0
 }
 ' tmp_schecked.tsv)

# Prepare input for liftover
# Extract relevant columns for target file, and create a BED file containing them. This will be the file fed to the liftover script 
 echo Preparing input for liftover...


 cat tmp_schecked.tsv | awk  -v snpidcol="$SNPIDCOL" -v chrcol="$CHRCOL" -v bpcol="$BPCOL" 'BEGIN{FS=OFS="\t"}{print "chr"$chrcol, $bpcol, (($bpcol + 1)), $snpidcol }' | tail -n+2 > ${FILEBASENAME}.bed

# In addition, we'll need to edit the file to make it "mergeable", since merge command will naturally try to merge using the first column in each file. We'll guarantee that the first column in both files to be SNPID.
 cp tmp_schecked.tsv tmp_formerging0.tsv

# We reorder the columns only if SNPID is not the first one
 if [[ "$SNPIDCOL" -gt 1 ]]; then
 	paste <(cut -f"${SNPIDCOL}" tmp_formerging0.tsv ) <(cut -f1-$((SNPIDCOL - 1)),$((SNPIDCOL + 1))- tmp_formerging0.tsv) > tmp_formerging1.tsv 
 	rm tmp_formerging0.tsv
 else
 	mv tmp_formerging0.tsv tmp_formerging1.tsv
 fi

# Check build
# echo Checking build...I will try to use SNPIDs to check for build. 
 echo Checking build...Using cupcake manifest SNPs to compare genomic coordinates. 
 	awk -v chrcol="$CHRCOL" -v bpcol="$BPCOL" 'BEGIN{FS="\t";OFS=":"}NR>1{ print $chrcol, $bpcol}' tmp_schecked.tsv > tmp_file_chrbp.tsv
 	matches18=$(grep -wf hg18_manifest.txt -c tmp_file_chrbp.tsv)
 	matches19=$(grep -wf hg19_manifest.txt -c tmp_file_chrbp.tsv)
 	matches38=$(grep -wf hg38_manifest.txt -c tmp_file_chrbp.tsv)
 	if [[ $matches38 -gt $matches19 ]] && [[ $matches38 -gt $matches18 ]];then
 		CHOSEN_BUILD=7
 	elif [[ $matches19 -gt $matches38 ]] && [[ $matches19 -gt $matches18 ]]; then
 		CHOSEN_BUILD=5
 	else 	CHOSEN_BUILD=3
 	fi
 	echo "$f had $matches18 matches with hg18, $matches19 with hg19, and $matches38 with hg38 for manifest SNPs."
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
 	echo "$f" seems to be in hg18.
 	echo Liftovering...
 	liftOver "${FILEBASENAME}".bed hg18ToHg38.over.chain.gz "${FILEBASENAME}"-lo-output.bed "${FILEBASENAME}"-unlifted.bed
 elif [ $CHOSEN_BUILD -eq 5 ]; then
 	echo "$f" seems to be in hg19.
 	echo Liftovering...
 	liftOver "${FILEBASENAME}".bed hg19ToHg38.over.chain.gz "${FILEBASENAME}"-lo-output.bed "${FILEBASENAME}"-unlifted.bed
 elif [ $CHOSEN_BUILD -eq 7 ]; then
 	echo "$f" is in hg38 already, skipping liftover step...
 	gzip < tmp_schecked.tsv  > "${FILEBASENAME}"-hg38.tsv.gz
 	continue # Remove/rethink this break when the whole pipeline is a single loop and there are more steps afterwards.
 else
 	echo Sorry, something wrong happened at the build selection stage and I could not identify the build for "$f".
 	echo "$f" could not be overlifted
 	continue 
 fi
 
 awk 'BEGIN{FS=OFS="\t"}{sub("chr", "",$1); print $4,$1,$2}' "${FILEBASENAME}"-lo-output.bed | sed '1i SNPID\tCHR38\tBP38' > "${FILEBASENAME}"-lo-output2.bed
 join -a1 -e'NA' -t $'\t' --nocheck-order -o auto "${FILEBASENAME}"-lo-output2.bed tmp_formerging1.tsv | gzip  > "${FILEBASENAME}"-hg38.tsv.gz
 echo "$f" suscessfully lifted over to hg38 build!
 
 rm  tmp_formerging1.tsv tempcolcheck.txt tmp_schecked.tsv *.bed



done
rm rs_manifest.txt hg18_manifest.txt hg19_manifest.txt hg38_manifest.txt

