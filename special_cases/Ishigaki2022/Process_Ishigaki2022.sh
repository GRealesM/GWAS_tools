# Processing Ishigaki 2022 RA for the basis

# Dw from paper repo
wget https://data.cyverse.org/dav-anon/iplant/home/kazuyoshiishigaki/ra_gwas/ra_gwas-10-28-2021.tar
tar -xf ra_gwas-10-28-2021.tar
mv 10-28-2021/* ./ # Bring files to current directory

# We'll focus on datasets with all individuals (although there are seropositive subsets)

Rscript --vanilla Process_Ishigaki2022.R EAS_all_auto-10-2021.txt.gz RA_Ishigaki_36333501_1-raw.tsv.gz
Rscript --vanilla Process_Ishigaki2022.R EUR_all_auto-10-2021.txt.gz RA_Ishigaki_36333501_2-raw.tsv.gz
Rscript --vanilla Process_Ishigaki2022.R Trans_all_auto-10-2021.txt.gz RA_Ishigaki_36333501_3-raw.tsv.gz

# Pipeline
/home/gr440/rds/rds-cew54-basis/GWAS_tools/01-Pipeline/pipeline_v5.3.2_beta.sh -f RA_Ishigaki_36333501_1-raw.tsv.gz
/home/gr440/rds/rds-cew54-basis/GWAS_tools/01-Pipeline/pipeline_v5.3.2_beta.sh -f RA_Ishigaki_36333501_2-raw.tsv.gz
/home/gr440/rds/rds-cew54-basis/GWAS_tools/01-Pipeline/pipeline_v5.3.2_beta.sh -f RA_Ishigaki_36333501_3-raw.tsv.gz

# Reduce for blood cell basis v2
Rscript --vanilla /home/gr440/rds/rds-cew54-basis/GWAS_tools/04-Reduction/Reducing_for_bloodcellbasis_v2.R RA_Ishigaki_36333501_1-hg38.tsv.gz
Rscript --vanilla /home/gr440/rds/rds-cew54-basis/GWAS_tools/04-Reduction/Reducing_for_bloodcellbasis_v2.R RA_Ishigaki_36333501_2-hg38.tsv.gz
Rscript --vanilla /home/gr440/rds/rds-cew54-basis/GWAS_tools/04-Reduction/Reducing_for_bloodcellbasis_v2.R RA_Ishigaki_36333501_3-hg38.tsv.gz

# Move to 02-Processed
mv *-hg38.tsv.gz ~/rds/rds-cew54-basis/02-Processed/