# Processing for Deakin et al. 2022 JDM dataset

# Several things that don't look right: no BETA (but SE), and no REF allele
sed -i -e '1s/chromosome/CHR/' -e '1s/variant_id/SNPID/' -e '1s/base_pair_location/BP/' -e '1s/effect_allele/ALT/' -e '1s/odds_ratio/OR/' -e '1s/standard_error/SE/' -e 's/p_value/P/' JDM_Deakin_35094092_1-raw.tsv
gzip JDM_Deakin_35094092_1-raw.tsv
# Check build (hg19!)
Rscript Deakin_2022.R 
# Also, not enough Basis coverage, but we'll project it anyway

zcat JDM_Deakin_35094092_1-raw.tsv.gz | sed -e '1s/CHR/CHR19/' -e '1s/BP/BP19/' > tmp.tsv; gzip tmp.tsv; mv tmp.tsv.gz JDM_Deakin_35094092_1-raw.tsv.gz

../GWAS_tools/00b-Missing_one_allele/fix-one-allele-problem.sh JDM_Deakin_35094092_1-raw.tsv.gz
zcat JDM_Deakin_35094092_1-fixed.tsv.gz | sed -e '1s/CHR19/CHR/' -e '1s/BP19/BP/' > tmp.tsv; gzip tmp.tsv; mv tmp.tsv.gz JDM_Deakin_35094092_1-fixed.tsv.gz
../GWAS_tools/01-Pipeline/pipeline_v5.3.2_beta.sh -f JDM_Deakin_35094092_1-fixed.tsv.gz
Rscript --vanilla ../GWAS_tools/04-Reduction/Reducing_for_IMDbasis.R JDM_Deakin_35094092_1-hg38.tsv.gz
Rscript --vanilla ../GWAS_tools/04-Reduction/Reducing_for_bloodcellbasis_v2.R JDM_Deakin_35094092_1-hg38.tsv.gz
mv JDM_Deakin_35094092_1-hg38.tsv.gz ../02-Processed/
Rscript ../GWAS_tools/05-Projection/Projecting_dbs_IMD_basis.R
Rscript ../GWAS_tools/05-Projection/Projecting_dbs_cell_basis_v3.R