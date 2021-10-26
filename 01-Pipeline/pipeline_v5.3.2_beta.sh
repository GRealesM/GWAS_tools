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
##

##################################################
## TARGET FILE REQUIREMENTS
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

##################################################
## REQUIRED FILES
##################################################

# This pipeline requires the following files to work
#
# - full_pipeline_vx.x.sh - This file
# - liftOver - UCSC liftOver executable
# - hg18ToHg38.over.chain.gz - hg18 to hg38 build Chain file
# - hg19ToHg38.over.chain.gz - hg19 to hg38 build Chain file
# - Manifest_build_translator.tsv - IMD basis manifest SNP genomic coordinates hg18, hg19, and hg38 builds
# - SE_Calc_pipeline.R - R script for calculating SE when absent from file
#

##################################################
## CHANGELOG from v4.1 
##################################################

# **v4.1** 
# 
# * I introduced a 's///' to remove carriage returns from files generated in Windows
# * Pipeline now joins liftovered files by keeping only SNPs that succesfully liftovered, as the previous behaviour (keeping all SNPs) was introducing way too many mismatches, including for SNPs that did in fact match. This may result in output files with less lines than input files.
# 
# **v4.4**
# 
# * Dictionary improved with more terms.
# * Pipeline now checks if BETA has more than 50% of SNPs as "NA", and recalculates from OR if it exists and has <50% NA.
# 
# **v4.5**
# 
# * I included a step in the final step (join, recompress, and save) to check and remove duplicated lines in the file. Because of the join command, some lines were systematically duplicated. This didn't have a huge effect on subsequent steps, but added innecessary lines.
# 
# **v4.6**
# 
# * SNPID is necessary for liftover step, but it's not absolutely required for it to contain the corresponding rs, so instead of throwing an error and jumping to next file when SNPID column is missing, it creates a new one from CHR and BP columns, with the CHR_BP format. 
# * Temporarily included A1 => REF and A2 => ALT to process Hoglund files.
# * Added new keys to the dictionary corresponding to Hoglund keys (eg. effB for BETA and se_effB for SE).
# 
# **v4.7**
# 
# * Fixed a bug that prevented cleaning-up of files when file is in hg38 build.
# 
# **v4.8**
# 
# * Included some terms in the dictionary to correctly process COVID-19 HGI datasets.
# 
# **v4.9**
# 
# * Included some terms in the dictionary to correctly process COVID-19 Rivas datasets.
# * Removed A1 and A2, included for specific purposes, from the dictionary.
# * Now once original build is identified, column names change to reflect it (eg. CHR/BP -> CHR19/BP19), rather than leave them untouched, as it did before.
# 
# **v4.10**
# 
# * New way to extract file base names, since now some of our datasets include dots in the trait name.
# 
# **v5.0**
# 
# * New post-liftover merging method, using awk instead of join. This result in seamless merge, keeping all original rows, using NA for missing hg38 coordinates, but without loss of rows due to faulty match, as it happened with the previous version.
# 
# **v5.1** 
#
# * Adapted to be call from any directory, and act over all .tsv.gz files in the directory from which it was called
#
# **v5.2**
#
# * Added a functionality to make allele columns uppercase
# * Added empty field to be recognised as missing data when checking the amount of missing data in key columns.
# * Added Help menu. 
# * Added "effect_allele_frequency" to the dictionary, to be replaced by ALT_FREQ.
# * Fixed bug that would make liftover merging fail (ie. printing all SNPs at the same position) when there's a SNPID column but it contains empty values.
#
# **v5.3**
#
# * Added the possibility to choose specific files, instead of blanket "all suitable files in directory"
# * Added creation of temporary directories, so pipeline can be run in parallel without fear to mix up temporary files.
# * Fixed a bug that would double columns in files with no SNPID column.
# * Changed the order of steps in liftOver step to save time when the file is in hg38 already.
#

shopt -s nullglob

version='5.3_beta'

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
	echo "GWAS Summary Statistics Harmonization Pipeline. Version $version"
	echo "This script is intended to harmonize GWAS summary statistics to a standard format, and includes the following steps:"
	echo "(1) - Detection of file separator and replacement by tabs."
	echo "(2) - Detection of column names and substitution by standard ones."
	echo "(3) - Sanity check of minimum columns, and computation of BETA when missing."
	echo "(4) - Detection of genome build and liftover from hg18 or hg19 to hg38."
	echo ""
	echo "This pipeline does not require any arguments, as it will work on all *tsv.gz files in the directory it was called. However, single files can now be specified."
	echo "The summary statistic datasets should have the following characteristics:"
	echo ""
	echo " - Compressed (in GZ format)."
	echo " - Fields separated by spaces, tabs, commas, or semicolons."
	echo " - Contains information on:"
	echo "	 * Chromosome (CHR)"
	echo "	 * Position (BP)"
	echo "	 * SNP id (SNPID)"
	echo "	 * Reference allele (REF)"
	echo "	 * Alternative, or effect allele (ALT)*"
	echo "	 * P-value for association (P)"
	echo "	 * Odds-ratio (OR) or BETA (BETA) "
	echo " 	 * Standard Error of the log OR (SE) "
	echo " Note: Tags in parentheses are suggestions, since the program will try to infer them from a number of possible tag names 	(Open file to see dictionary)."
	echo " * We've seen wild heterogeneity in Effect/Reference allele nomenclature (eg. A1/A2 meaning either REF/ALT or ALT/REF, depending on the study)"
	echo " so I recommend to explicitly rename them prior to running the pipeline."
	echo " - Must be in either GRCh36/hg18, GRCh37/hg19, or GRCh38/hg38 builds."
	echo ""
   	echo "options:"
        echo "h     Print this Help."
   	echo "V     Print software version and exit."
	echo "f     Optional: Specify single file instead of all .tsv.gz in the directory"
}

Version(){
	echo "v$version"
}

# Get the options
while getopts ":f:Vh" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      f) # Specifies .tsv.gz file, if desired
	 array=${OPTARG}
	 ;;
      V) # display version
	 Version
	 exit;;
      \?) # incorrect option
         echo "Error: Invalid option"
	 Help
         exit;;
   esac
done

echo "GWAS summary statistics processing pipeline. Version $version."

# Check if argument was provided, and choose all (classic option) otherwise.
if [ -z "$array" ]; then
	echo "Running pipeline for all .tsv.gz files in the directory"
	array=(*.tsv.gz) # Initial file list, make sure to have all in tsv.gz format, since chain files are .gz, so including only .gz will try to process those, too :/
	else
	echo "Running pipeline for $array."
fi

#############################################
## GLOBAL VARIABLES (Generated once per run)
#############################################

scriptpath=(~/rds/rds-cew54-basis/GWAS_tools/01-Pipeline/)



#######################
## MAIN LOOP
#######################



for f in "${array[@]}";
do

 FILEBASENAME=$(echo "$f" | sed -e 's/.tsv.gz//' -e 's/-.*//'); # Take as file base name the first dash or file extension.

## Create a temporary directory, copy the target file to it, cd to it, and then run the
# rest of the processing in it
# Create a temporary directory and store its name in a variable 
TEMPLATE_PREFIX='tmp' # prefix of your new tempdir template
TEMPLATE_RANDOM='XXXX' # Increase the Xs for more random characters
TEMPLATE=${TEMPLATE_PREFIX}.${TEMPLATE_RANDOM}

# create the tempdir using your custom $TEMPLATE, which may include
# a path such as a parent dir, and assign the new path to a var
TMPDIR=$(mktemp -d $TEMPLATE)

# Bail out if the temp directory wasn't created successfully.
if [ ! -e $TMPDIR ]; then
    >&2 echo "Failed to create temp directory"
    exit 1
fi

# Make sure it gets removed even if the script exits abnormally.
#trap "exit 1"           HUP INT PIPE QUIT TERM
#trap 'rm -rf ../"$TMPDIR"' EXIT

cp "$f" "$TMPDIR"
cd "$TMPDIR"

## IMD Basis SNP ids and coordinates for the 3 available builds, for build detection

awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1}' "$scriptpath"Manifest_build_translator.tsv > rs_manifest.txt
awk 'BEGIN{FS="\t";OFS=":"} NR>1 {print $2,$3}' "$scriptpath"Manifest_build_translator.tsv > hg18_manifest.txt
awk 'BEGIN{FS="\t";OFS=":"} NR>1 {print $4,$5}' "$scriptpath"Manifest_build_translator.tsv > hg19_manifest.txt
awk 'BEGIN{FS="\t";OFS=":"} NR>1 {print $6,$7}' "$scriptpath"Manifest_build_translator.tsv > hg38_manifest.txt

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


### Only for the KUBO files processing, I temporarily added Allele1 = REF and Allele2 = ALT. This may not be the case for other files, so remember to remove them later!

 echo Reheading...
 	head -n1 tmp_fs_file.tsv | awk -F'\t' '{
 	gsub(/\<Chr\>|\<chromosome\>|\<Chromosome\>|\<chr\>|\<Chr_ID\>|\<hg18chr\>|\<CHROMOSOME\>|#chrom\>|#CHROM\>|\<chrom\>|#CHR\>/,"CHR");
 	gsub(/\<Pos\>|\<base_pair_location\>|\<BP\>|\<BP\(hg19\)\>|\<Position\>|\<POS\>|\<pos\>|\<Chr_Position\>|\<bp\>|\<position\>|\<Position\(hg19\)\>|\<POSITION\>|\<bp_hg19\>|\<Coordinate\>|\<chrloc\>/,"BP");
 	gsub(/\<íd\>|\<id\>|\<ID\>|\<variant_id\>|\<MarkerName\>|\<SNP\>|\<rsid\>|\<rsids\>|\<SNP_Name\>|\<snp\>|\<snpid\>|\<SNP_ID\>|\<rsID\>|#SNPID\>|\<rs_number\>|\<RSID\>|\<rs\>|\<db_SNP_RS_ID\/Marker\>|\<dbSNP_RS_ID\>|\<Variant\>/,"SNPID");
 	gsub(/\<OtherAllele\>|\<reference_allele\>|\<Ref_Allele\>|\<OTHER_ALLELE\>|\<other_allele\>|\<A2_other\>|\<NEA\>|\<Ref_Allele\>|\<Ref\>|\<ref\>|\<Allele1\>/,"REF");
 	gsub(/\<effect_allele\>|\<Effect_Allele\>|\<EffectAllele\>|\<A1_effect\>|\<RISK_ALLELE\>|\<EA\>|\<Risk_Allele\>|\<EFFECT_ALLELE\>|\<Alt\>|\<alt\>|\<Allele2\>/,"ALT");
 	gsub(/\<Beta\>|\<beta\>|\<Effect\>|\<effect\>|\<EFFECT\>|\<beta_SNP_add\>|\<EFFECT_ALT\>|\<effB\>|\<all_inv_var_meta_beta\>/,"BETA");
 	gsub(/\<standard_error\>|\<StdErr\>|\<stderr\>|\<sebeta_SNP_add\>|\<se\>|\<STDERR\>|\<sebeta\>|\<se_effB\>|\<all_inv_var_meta_sebeta\>|\<LOG\(OR\)_SE\>/,"SE");
 	gsub(/\<odds_ratio\>|\<Odds_ratio\>|\<or\>|\<OddsRatio\>|\<OR\(A1\)\>|\<ORX\>/,"OR");
 	gsub(/\<p_value\>|\<P.value\>|\<pvalue\>|\<P-value\>|\<pval\>|\<p.value\>|\<Pval\>|\<PVALUE\>|\<Pvalue\>|\<P_VALUE\>|\<P-val\>|\<p\>|\<All.p.value\>|\<P_value\>|\<p-value\>|\<GC-adjusted_P_\>|\<Chi-Squared__P\>|\<P1df\>|\<all_inv_var_meta_p\>/,"P");
 	gsub(/\<Log10p\>/,"LOG10P");
 	gsub(/\<_-log10_p-value\>/,"-LOG10P");
 	gsub(/\<effect_allele_frequency\>|<\maf\>|<\MAF\>/,"ALT_FREQ");
 # Caution! Sometimes "other_allele" means effect allele, check papers prior to run the script, and pre-rename accordingly.
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
 	gsub(/\<Rsq\>/,"RSQ");
 	gsub(/\<MARKER\>|\<íd\>|\<Chr\:Position\>/,"CHR:BP");
	gsub(/\<Zscore\>|\<ZSCORE\>|\<Z_STAT\>/,"Z");print}' > tmp_reheaded_file.tsv
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
  	PCOL=$(awk -F'\t' '{for(i=1;i<=NF;i++) {if($i == "P") printf(i)	} exit 0}' tmp_reheaded_file.tsv)
	NAP=$(echo "scale=2;$(awk -v pcol="$PCOL" 'BEGIN{FS=OFS="\t"}NR>1{print $pcol}' tmp_reheaded_file.tsv | grep -cP 'NA|^$') / ($(cat tmp_reheaded_file.tsv | wc -l)-1)" | bc)
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
 	BETACOL=$(awk -F'\t' ' { for(i=1;i<=NF;i++) { if($i == "BETA") printf(i) } exit 0 }' tmp_reheaded_file.tsv)
	NABETA=$(echo "scale=2;$(awk -v betacol="$BETACOL" 'BEGIN{FS=OFS="\t"}NR>1{print $betacol}' tmp_reheaded_file.tsv | grep -cP 'NA|^$') / ($(cat tmp_reheaded_file.tsv | wc -l)-1)" | bc)
	if [[ 1 -eq $(echo  "$NABETA > $maxNA" | bc -l) ]]; then
		echo BETA column seems to have more NAs than what is acceptable. Please check your file.
		echo Continuing with next file...
		continue
	fi	
fi

# End of first stage. Up to here we have only checked if the very minimum columns required are present. Some of the columns (ie BETA and SE) can be calculated from available data

# FIXING ALLELE CASE
# Some files come with alleles in lower case, which is different from our standard, let's fix that
echo Ensuring allele columns are uppercase...
REFCOL=$(awk -F'\t' '{for(i=1;i<=NF;i++) { if($i == "REF") printf(i) } exit 0}' tmp_reheaded_file.tsv)
ALTCOL=$(awk -F'\t' '{for(i=1;i<=NF;i++) { if($i == "ALT") printf(i) } exit 0}' tmp_reheaded_file.tsv)
awk -v refcol="$REFCOL" -v altcol="$ALTCOL" 'BEGIN{FS=OFS="\t"}NR>1{$refcol=toupper($refcol); $altcol=toupper($altcol)}1' tmp_reheaded_file.tsv > tmp_schecked.tsv

rm tmp_reheaded_file.tsv

# We can calculate BETA from OR (if OR column exists and is not all NA).
if [[  "$minOR" == 1 && "$minBETA" == 0 ]]; then	
 	ORCOL=$(awk -F'\t' ' { for(i=1;i<=NF;i++) {if($i == "OR") printf(i)} exit 0}' tmp_schecked.tsv)
	NAOR=$(echo "scale=2;$(awk -v orcol="$ORCOL" 'BEGIN{FS=OFS="\t"}NR>1{print $orcol}' tmp_schecked.tsv | grep -cP 'NA|^$') / ($(cat tmp_schecked.tsv | wc -l)-1)" | bc)
	if  [ 1 -eq $(echo  "$NAOR > $maxNA" | bc -l) ]; then
		echo OR column seems to have more NAs than what is acceptable. Please check your file.
		echo Continuing with next file...
		continue	
	else 	
		echo "$f" lacks BETA but has OR. BETA will be calculated 
 		awk -v orcol="$ORCOL" 'BEGIN{FS=OFS="\t"} {print $0,log($orcol)}' tmp_schecked.tsv | sed -e 's///' -e '1s/-inf/BETA/' > tmp.tsv && mv tmp.tsv tmp_schecked.tsv
	fi
fi


# If BETA column exist but has too many NAs, and OR exist, we can recalculate it
if [[ "$minBETA" == 1 && 1 -eq $(echo  "$NABETA > $maxNA" | bc -l) ]]; then
	if [[ "$minOR" == 0 ]]; then
		echo BETA column seems to have more NAs than what is acceptable, and there is no OR column to recalculate it from. Please check your file.
		echo Continuing with next file...
		continue
	elif [[ "$minOR" == 1 ]]; then
		NAOR=$(echo "scale=2;$(awk -v orcol="$ORCOL" 'BEGIN{FS=OFS="\t"}NR>1{print $orcol}' tmp_schecked.tsv | grep -cP 'NA|^$') / ($(cat tmp_schecked.tsv | wc -l)-1)" | bc)
		if [[ 1 -eq $(echo  "$NAOR > $maxNA" | bc -l) ]]; then
			echo BETA and OR both have more NAs that is acceptable. Please check your file.
			echo Continuing with next file...
			continue
		else
			echo BETA column has more NA than acceptable, but file has OR. BETA will be calculated.
			cut -f"$BETACOL" --complement tmp_schecked.tsv | awk -v orcol="$ORCOL" 'BEGIN{FS=OFS="\t"} {print $0,log($orcol)}' | sed -e 's///' -e '1s/-inf/BETA/' > tmp.tsv && mv tmp.tsv tmp_schecked.tsv
		fi
	fi
fi


# We can calculate SE if we have BETA and Z or P values
if [[ "$minSE" == 0 ]]; then
	echo "$f" seems to lack SE. It will be calculated using BETA and P values. This will calculate Z scores as a side effect too.
	Rscript "$scriptpath"SE_Calc_pipeline.R
	mv tmp.tsv tmp_schecked.tsv
fi


# A couple annoying things that we should take care of are
# (1) In some files the X chromosome will be encoded as 23. Apparently the almighty liftOver executable is not very happy with this notation
# (2) Likewise, in some files, BP are expressed as scientific notation, which also makes liftOver cringe in disgust.
# (3) In addition, in some files CHR is presented as chrX, which will likely upset liftOver too. We need to remove it.
# so we'll need to add a little workaround
# (4) Lastly, we'll check if the SNPID column exists, if not, we'll create one from CHR and BP. For it not to enter in conflict with existing CHR:BP columns, we'll use an underscore instead of a semicolon.
CHRCOL=$(awk -F'\t' '{ for(i=1;i<=NF;i++) { if($i == "CHR") printf(i) } exit 0 }' tmp_schecked.tsv)
BPCOL=$(awk -F'\t' ' { for(i=1;i<=NF;i++) { if($i == "BP")  printf(i) } exit 0 }' tmp_schecked.tsv)
REFCOL=$(awk -F'\t' '{ for(i=1;i<=NF;i++) { if($i == "REF") printf(i) } exit 0 }' tmp_schecked.tsv)
ALTCOL=$(awk -F'\t' ' { for(i=1;i<=NF;i++) { if($i == "ALT")  printf(i) } exit 0 }' tmp_schecked.tsv)


awk -v chrcol="$CHRCOL" -v bpcol="$BPCOL" 'BEGIN{FS=OFS="\t"}{sub(/23/,"X",$chrcol); sub("chr","",$chrcol)}NR>1{$bpcol=sprintf("%i", $bpcol)}1' tmp_schecked.tsv > tmp.tsv && mv tmp.tsv tmp_schecked.tsv


if [[ "$minSNPID" == 0 ]]; then
 	echo "$f" seem to lack SNPIDs. This is not ideal, but I will create one using CHR and BP
 	awk -v chrcol="$CHRCOL" -v bpcol="$BPCOL" -v refcol="$REFCOL" -v altcol="$ALTCOL" 'BEGIN{FS=OFS="\t"} {print $0,$chrcol"_"$bpcol"_"$refcol"_"$altcol}' tmp_schecked.tsv | sed -e 's///' -e '1s/CHR_BP_REF_ALT/SNPID/' > tmp.tsv && mv tmp.tsv tmp_schecked.tsv

fi

 SNPIDCOL=$(awk -F'\t' '{for(i=1;i<=NF;i++) {if($i == "SNPID") printf(i)} exit 0 }' tmp_schecked.tsv)
 NASNPID=$(awk -v snpidcol="$SNPIDCOL" 'BEGIN{FS=OFS="\t"}NR>1{print $snpidcol}' tmp_schecked.tsv | grep -cP 'NA|^$')

if [[ "$NASNPID" != 0 ]]; then
	echo "SNPID column seems to have some missing data. We'll use CHR:BP format instead"
	awk -v snpidcol="$SNPIDCOL" -v chrcol="$CHRCOL" -v bpcol="$BPCOL" -v refcol="$REFCOL" -v altcol="$ALTCOL" 'BEGIN{FS=OFS="\t"}($snpidcol == "NA" || $snpidcol == "") {$snpidcol=$chrcol"_"$bpcol"_"$refcol"_"$altcol}1' tmp_schecked.tsv > tmp.tsv && mv tmp.tsv tmp_schecked.tsv
        NASNPID=$(awk -v snpidcol="$SNPIDCOL" 'BEGIN{FS=OFS="\t"}NR>1{print $snpidcol}' tmp_schecked.tsv | grep -cP 'NA|^$')
	if [[ "$NASNPID" != 0 ]]; then
		echo "SNPID column seems to still have some missing values, removing them"
		awk -v snpidcol="$SNPIDCOL" 'BEGIN{FS=OFS="\t"} $snpidcol !~ /NA|^$/ ' tmp_schecked.tsv > tmp.tsv && mv tmp.tsv tmp_schecked.tsv 
	fi
fi

echo Column sanity check OK. 
	

###################
## LIFTOVER STAGE
###################


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
	sed -e '1s/\<CHR\>/CHR18/' -e '1s/\<BP\>/BP18/' tmp_formerging1.tsv > tmp.tsv && mv tmp.tsv tmp_formerging1.tsv
 	echo Liftovering...
 	"$scriptpath"liftOver "${FILEBASENAME}".bed "$scriptpath"hg18ToHg38.over.chain.gz "${FILEBASENAME}"-lo-output.bed "${FILEBASENAME}"-unlifted.bed

elif [ $CHOSEN_BUILD -eq 5 ]; then
 	echo "$f" seems to be in hg19.
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
	sed -e '1s/\<CHR\>/CHR19/' -e '1s/\<BP\>/BP19/' tmp_formerging1.tsv > tmp.tsv && mv tmp.tsv tmp_formerging1.tsv
 	echo Liftovering...
 	"$scriptpath"liftOver "${FILEBASENAME}".bed "$scriptpath"hg19ToHg38.over.chain.gz "${FILEBASENAME}"-lo-output.bed "${FILEBASENAME}"-unlifted.bed

elif [ $CHOSEN_BUILD -eq 7 ]; then
 	echo "$f" is in hg38 already, skipping liftover step...
 	sed -e '1s/\<CHR\>/CHR38/' -e '1s/\<BP\>/BP38/' tmp_schecked.tsv | gzip  > ../"${FILEBASENAME}"-hg38.tsv.gz
 	cd .. && rm -rf "$TMPDIR"
 	continue # Remove/rethink this break when the whole pipeline is a single loop and there are more steps afterwards.
 else
 	echo Sorry, something wrong happened at the build selection stage and I could not identify the build for "$f".
 	echo "$f" could not be overlifted
 	continue 
 fi
 
 
 awk 'BEGIN{FS=OFS="\t"}{sub("chr", "",$1); print $4,$1,$2}' "${FILEBASENAME}"-lo-output.bed | sed '1i SNPID\tCHR38\tBP38' > "${FILEBASENAME}"-lo-output2.bed
 # New joining command from 5.0 on
 awk -v OFS='\t' 'NR==FNR{a1[$1]=$2; a2[$1]=$3;next};{ if ($1 in a1) print $0, a1[$1], a2[$1]; else print $0, "NA","NA"}' "${FILEBASENAME}"-lo-output2.bed tmp_formerging1.tsv | gzip  > ../"${FILEBASENAME}"-hg38.tsv.gz
#  snpsbeforeliftover=$(echo ""$(cat tmp_formerging1.tsv | wc -l)" - 1" | bc)
#  snpsafterliftover=$(echo ""$(zcat "${FILEBASENAME}"-hg38.tsv.gz | wc -l)" -1" | bc)
#  snpsdifference=$(echo "$snpsbeforeliftover" - "$snpsafterliftover" | bc)
#  echo ""$f" had "$snpsbeforeliftover" SNPs. After liftover it now has "$snpsafterliftover" ("$snpsdifference" less)."
 echo "$f" suscessfully lifted over to hg38 build!
 
# rm  tmp_formerging1.tsv tempcolcheck.txt tmp_schecked.tsv *.bed rs_manifest.txt hg18_manifest.txt hg19_manifest.txt hg38_manifest.txt

 cd .. && rm -rf "$TMPDIR"

done

