# Fixing PBC_Nakamura_23000144_1
# 2020-03-05
# I downloaded the file again, and it has some particularities regarding alleles. BETA/OR refers to the minor allele, which changes depending on the SNP, there's a column specifying which one is the minor allele. Let's fix it!
# The original file can be downloaded here: ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/NakamuraM_23000144_GCST001685/hum0076_1stgwas_160916.csv
# The file in 97-Archive is already processed by the code below, but not pipelined.

library(data.table)

dt <- fread("PBC_Nakamura_23000144_1-raw.tsv")
# Remove [ and ] 
dt[, `Reference Alleles A/B`:=gsub("\\[|\\]", "", `Reference Alleles A/B`)]
# Separate alleles into two columns
dt[, c("REF", "ALT") := tstrsplit(`Reference Alleles A/B`, "\\/")]
# calculate BETA from OR
dt[, BETA := log(`Odds Ratio (Minor Allele)`)]
# Invert beta sign depending on which allele it referst to. Since it refers to minor allele, and we took B to be our ALT, we'll need to switch signs only when Minor allele = A.
dt[`Minor Allele` == "A", `:=` (BETA=BETA*-1, `Minor Allele Frequency`= 1-`Minor Allele Frequency`)]
# Remove unneeded columns
dt[, c("Reference Alleles A/B", "Odds Ratio (Minor Allele)", "OR Lower Confidence Bound (Minor)", "OR Upper Confidence Bound (Minor)", "Minor Allele") := NULL]
setnames(dt, "Minor Allele Frequency", "ALT_Freq")
fwrite(dt, "PBC_Nakamura_23000144_1-raw.tsv.gz", sep = "\t", na = "NA", row.names = FALSE, quote = FALSE)
