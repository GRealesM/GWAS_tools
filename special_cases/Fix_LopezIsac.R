# Fixing SSC_LopezIsac
# This file has only one allele reported. I took the SNP alleles from 1000genomes to take the ther one

library(data.table)

lopezisac <- fread("SSC_LopezIsac_31672989_1-raw.tsv.gz")
alleles <- fread("tomerge_1kg_LopezIsac_alleles.txt", sep = ":")
lopezisac <- merge.data.table(lopezisac, alleles[,c(1, 3:4)], by.x = "SNP", by.y = "V1")
lopezisac$A2 <- ifelse(lopezisac$A1 == lopezisac$V3, lopezisac$V4, lopezisac$V3)
lopezisac[,c("V3", "V4") := NULL]
lopezisac <- unique(lopezisac)
fwrite(lopezisac, "SSC_LopezIsac_31672989_1-fixed.tsv.gz", sep = "\t", row.names = FALSE, quote = FALSE)
