# Selecting and downloading Orrú et al., 2020

# Guillermo Reales
# 2023/03/17

# Load packages
library(data.table)
setDTthreads(15)


f <- fread("https://www.ebi.ac.uk/gwas//api/search/summaryStatistics/download")
o <- f[`PubMed ID` == 32929287]
o <- o[grepl("count", `Reported Trait`, ignore.case = TRUE)][ grepl("lymphocyte", `Trait(s)`)]
s <- o[c(5:8, 25:26, 38:40, 42:43, 50, 62)] # Selected subset
s
#  1:                      Switched memory B cell Absolute Count
#  2:                    Plasma Blast-Plasma Cell Absolute Count
#  3:                               Memory B cell Absolute Count
#  4:                         Naive-mature B cell Absolute Count
#  5:                           Naive CD8+ T cell Absolute Count
#  6:                 Effector Memory CD8+ T cell Absolute Count
#  7:                  Central Memory CD4+ T cell Absolute Count
#  8:                           Naive CD4+ T cell Absolute Count
#  9:                 Effector Memory CD4+ T cell Absolute Count
# 10:                  Central Memory CD8+ T cell Absolute Count
# 11:               Resting CD4 regulatory T cell Absolute Count
# 12: Activated & secreting CD4 regulatory T cell Absolute Count
# 13:                    Unswitched memory B cell Absolute Count

# We'll download the GWAS catalog harmonised version of the files. Despite what it says in the "Trait URI" column, files all have EFO_0007937 in their URL
s[, URL:=paste0("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90001001-GCST90002000/", `Study Accession`, "/harmonised/", `PubMed ID`,"-",`Study Accession`, "-EFO_0007937.h.tsv.gz")]
nms  <- c("smBC","pbPC","mBC","nmBC","nCD8","emCD8","cmCD4","nCD4","emCD4","cmCD8","rCD4reg","asCD4reg","umBC")
s[, tname:=nms]
s[, fname:=paste0(nms, "_Orru_32929287_1-raw.tsv.gz")]

for(i in 1:13){
    message("Downloading ", s$tname[i], "...")
    system(paste0("wget ", s$URL[i], " -O ", s$fname[i]))
}


# Now we can prepare them for GWAS preparation pipeline

for(i in 1:13){
    message("Fixing headers in ", paste0(nms[i], "_Orru_32929287_1-raw.tsv.gz"))
    fs <- fread(paste0(nms[i], "_Orru_32929287_1-raw.tsv.gz"))
    setnames(fs, c("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele", "n", "effect_allele_frequency", "beta", "standard_error", "p_value", "odds_ratio"),
                 c("SNPID","CHR", "BP", "ALT", "REF", "N", "ALT_FREQ", "BETA", "SE", "P", "OR"))
    fwrite(fs, file = paste0(nms[i], "_Orru_32929287_1-hc.tsv.gz"), sep="\t", na = NA, quote=FALSE)
}

# Now apply pipeline

for(i in 2:13){
    message("Pipelining ", paste0(nms[i], "_Orru_32929287_1-hc.tsv.gz"))
    system(paste0("~/rds/rds-cew54-basis/GWAS_tools/01-Pipeline/pipeline_v5.3.2_beta.sh -f ", nms[i], "_Orru_32929287_1-hc.tsv.gz"))
}

# The paper states that all measurements were normalised, so no need to adjust beta and SE.
# Now we just need to reduce them (sbatch --array 1-13 slurm_redcb2_massive) and project them