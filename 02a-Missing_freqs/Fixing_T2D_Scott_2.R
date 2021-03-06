# Fixing T2D_Scott_
# Created 2020-03-15
# This file had the issue of not having SNPIDs. Since we don't *need* ids to process the file, we'll separate the CHR and POS, and keep the CHRPOS column as SNPID

library(data.table)

t2d <- fread("../../97-Archive/T2D_Scott_28566273_1.tsv.gz")
t2d[, c("CHR", "BP") := tstrsplit(`Chr:Position`, ":")]
t2d$`Chr:Position`  <- gsub(pattern=":", replacement="_", x=t2d$`Chr:Position`)
setnames(t2d, old = c("Chr:Position", "Allele1", "Allele2", "Effect", "StdErr"), new = c("SNPID", "ALT", "REF", "BETA", "SE"))
fwrite(t2d, "../../01-Renamed/T2D_Scott_28566273_1-fixed.tsv.gz", sep = "\t", row.names = FALSE, quote = FALSE)
