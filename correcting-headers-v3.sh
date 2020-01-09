#!/bin/bash

shopt -s nullglob
array=(*.gz)

for f in ${array[@]};
do
    FILEBASENAME=$(echo $f | cut -d- -f1); # Take as file base name everything before the first dash
    echo "Working on $FILEBASENAME file"
    zcat $f | awk -F"\t" '{\
	gsub(/\<Log10p\>/,"LOG10P");\
	gsub(/\<_-log10_p-value\>/,"-LOG10P");\
	gsub(/\<effect_allele\>|\<EffectAllele\>|\<A1_effect\>|\<RISK_ALLELE\>|\<EA\>|\<Risk_Allele\>|\<EFFECT_ALLELE\>|\<Alt\>/,"ALT");\
	gsub(/\<OtherAllele\>|\<reference_allele\>|\<OTHER_ALLELE\>|\<other_allele\>|\<A2_other\>|\<NEA\>|\<Ref_Allele\>|\<Ref\>/,"REF");\
# Caution! Sometimes "other_allele" means effect allele, check papers prior to run the script, and pre-rename accordingly.
	gsub(/\<Beta\>|\<beta\>|\<Effect\>|\<effect\>|\<EFFECT\>|\<sebeta_SNP_add\>|\<EFFECT_ALT\>/,"BETA");\
	gsub(/\<Pos\>|\<base_pair_location\>|\<BP\>|\<BP\(hg19\)\>|\<Position\>|\<POS\>|\<pos\>|\<Chr_Position\>|\<bp\>|\<position\>|\<Position\(hg19\)\>|\<POSITION\>|\<bp_hg19\>|\<Coordinate\>|\<chrloc\>/,"BP");\
	gsub(/\<Chr\>|chromosome\>|\<Chromosome\>|\<chr\>|\<Chr_ID\>|\<hg18chr\>|CHROMOSOME\>/,"CHR");\
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
	gsub(/\<id\>|\<variant_id|>|\<MarkerName\>|\<SNP\>|\<rsid\>|\<SNP_Name\>|\<snp\>|\<snpid\>|\<rsID\>|\<\\#SNPID\>|\<rs_number\>|\<RSID\>|\<rs\>|\<db_SNP_RS_ID\/Marker\>|\<dbSNP_RS_ID\>/,"SNPID");\
	gsub(/\<MARKER\>|\<íd\>|\<Chr\:Position\>/,"CHR:BP"); print}' |\
    gzip > ../03-Header_corrected/"${FILEBASENAME}-hc.tsv.gz"
done
echo "Done, your files are now at ... directory"
