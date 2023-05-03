# Process PAPS
setwd("/home/gr440/rds/rds-cew54-basis/01-Workshop")
library(data.table)
library(magrittr)

d <- fread("All_studies_EA_IT_NE_SP1_SP2_rsq9_RefAl.meta")
d[, A2:= sapply(strsplit(SNP,":"),`[`, 4)] # Extract A2 info from SNP column
names(d)[3:5] <- c("SNPID", "ALT", "REF")
setcolorder(d, c("CHR", "BP", "SNPID","REF", "ALT"))

fwrite(d, "PAPS_Casares_up_1-raw.tsv.gz", sep = "\t", na=NA, quote=F, nThread = 20)
