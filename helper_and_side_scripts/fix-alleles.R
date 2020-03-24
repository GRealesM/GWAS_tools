# Fixing files with only one allele reported
# I took the SNP alleles from 1000genomes to take the ther one

library(data.table)

alleles <- fread("Alleles.txt", sep = ":")

for(i in dir(pattern = "*.tsv.gz")){
  filetofix <- fread(i)
  filetofix <- merge.data.table(filetofix, alleles[,c(1, 3:4)], by.x = "SNP", by.y = "V1")
  filetofix$A2 <- ifelse(filetofix$A1 == filetofix$V3, filetofix$V4, filetofix$V3)
  filetofix[,c("V3", "V4") := NULL]
  filetofix <- unique(filetofix)
  fwrite(filetofix, i, sep = "\t", row.names = FALSE, quote = FALSE)
}
