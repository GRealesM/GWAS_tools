# Fixing files with only one allele reported
# I took the SNP alleles from 1000genomes to take the ther one

library(data.table)

alleles <- fread("alleles.txt", sep = ":", fill=TRUE)
alleles$coord  <- paste(alleles$V1, alleles$V2, sep=":")

files  <- dir(pattern="*tsv.gz")

for (file in files){

	filebasename <- strsplit(file, split=".", fixed = TRUE)[[1]][1]
	filetofix <- fread(file)
	message(file, " has ", nrow(filetofix), " SNPs.")
	filetofix$coord  <- paste(filetofix$CHR, filetofix$BP, sep=":")
	filetofix <- merge.data.table(filetofix, alleles[,c(6, 4:5)], by="coord", all.y=TRUE)
	# NOTE: Also, I renamed A1 to ALT, because that was the one-allele denomination in Traglia. Change it to suit your needs.
	setnames(filetofix, old="A1", new="ALT")
	filetofix$REF <- ifelse(filetofix$ALT == filetofix$V5, filetofix$V4, filetofix$V5)
	filetofix[,c("V4", "V5", "coord") := NULL]
	filetofix <- unique(filetofix)
	setcolorder(filetofix, c("SNP","CHR","BP", "REF","ALT", "BETA", "SE", "P"))
	message("Processed ", paste(filebasename, "-fixed.tsv.gz", sep="")," has ", nrow(filetofix), " SNPs. Writing...")
	fwrite(filetofix, paste(filebasename, "-fixed.tsv.gz", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
	message("Done!")

}
