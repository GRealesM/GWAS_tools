Kazuyoshi Ishigaki, MD. PhD
10-28-2021


I, GWAS summary statistics
[GWAS type]
1, Autosomal, all RA cases (seropositive + seronegative) vs controls
- Trans_all_auto-10-2021.txt.gz: Trans-ancestry meta-analysis (37 cohorts: 35,871 cases and 240,149 controls)
- EUR_all_auto-10-2021.txt.gz: EUR meta-analysis, autosomal (25 cohorts: 22,350 cases and 74,823 controls)
- EAS_all_auto-10-2021.txt.gz: EAS meta-analysis, autosomal (8 cohorts: 11,025 cases and 162,608 controls)

2, Autosomal, seropositive RA cases vs controls
- Trans_seroposi_auto-10-2021.txt.gz (37 cohorts: 27,448 cases vs 240,149 controls)
- EUR_seroposi_auto-10-2021.txt.gz (25 cohorts: 17,221 cases and 74,823 controls)
- EAS_seroposi_auto-10-2021.txt.gz (8 cohorts: 8,340 cases and 162,608 controls)

3, Chr X, all RA cases (seropositive + seronegative) vs controls
- Trans_all_chrx-10-2021.txt.gz

[File format (header)]
- SNP: chr_pos_ref_alt (hg19)
- Beta: effect size estimates (the effect allele is the alternative allele)
- SE: S.E. of the effect size estimate
- Pval: P value



II, PRS model
- prs_model_transgwas_impact.txt.gz: the PRS model which can be readily used for Plink2. This is the model with the best performance in all ancestries in our study (trans-ancestry GWAS + CD4T-Tbet-top5%)
- Important: this PRS model is based on the autosomal variants excluding the MHC region. Including the MHC region would surely improve the PRS performance.

[File format (header)]
- SNP: chr_pos_ref_alt (hg19)
- Effect_allele: effect allele
- Beta: effect size estimate


