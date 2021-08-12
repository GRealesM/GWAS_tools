# Inshaw dataset is a bit tricky, because it is a meta-analysis, it doesn't have clear REF/ALT distribution

library(data.table)

ds <- fread("~/rds/rds-cew54-wallace-share/Data/GWAS-summary/inshaw_t1d_1kg_imp.txt.gz")
names(ds)

table(ds$minor_allele == ds$alleleA)
# 
#   FALSE    TRUE 
# 7139088 1851095 
table(ds$minor_allele == ds$alleleB)
# 
#   FALSE    TRUE 
# 1851095 7139088 
# In most occasions, AlleleB matches the minor_allele. If we assume that the effect allele is usually the minor allele, we can say alleleA = REF and alleleB = ALT
# We choose CHR:BP:REF:ALT as SNPID since rsids has missing values, classified as "." which pipeline 5.1 won't deal well with 
setnames(ds, c("id","chromosome", "position","seall","beta","pmeta","alleleA", "alleleB", "ref"), c("SNPID","CHR", "BP", "SE", "BETA", "P", "REF", "ALT", "reference.source"))
# Remove duplicated SNPID too
ds <- unique(ds, by="SNPID")
fwrite(ds, "T1D_Inshaw_up_1-raw.tsv.gz", sep="\t")
