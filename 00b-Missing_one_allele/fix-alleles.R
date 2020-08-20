# Fixing files with only one allele reported
# I took the SNP alleles from 1000genomes to take the ther one

library(data.table)

switches <- c(A = "T", T = "A", C = "G", G = "C")

allele.diag <- function(x){

	flip.class <- ifelse(x$ALT == x$V5, 
			     "nochange",
			     ifelse(x$ALT == x$V4, 
				    "reverse",
				    ifelse(switches[x$ALT]== x$V5, 
					   "comp", 
					   ifelse(switches[x$ALT] == x$V4,
						  "revcomp",
						  "impossible"))))
	return(flip.class)
}

file <- commandArgs(trailingOnly = TRUE)
if(length(file) != 1) stop("An argument must be provided!")
alleles <- fread("alleles.txt", sep = ":", fill=TRUE)
alleles$coord  <- paste(alleles$V1, alleles$V2, sep=":")

filebasename <- strsplit(file, split="-|.tsv.gz")[[1]][1]
filetofix <- fread(file)
message(file, " has ", nrow(filetofix), " SNPs.")
filetofix$coord  <- paste(filetofix$CHR19, filetofix$BP19, sep=":")
filetofix <- merge.data.table(filetofix, alleles[,c(6, 4:5)], by="coord")
# Sometimes there are trialellic and tetra-allelic SNPs. Those are hard to make a decision on, so we'll remove them
filetofix <- filetofix[!duplicated(filetofix$coord),]
flip.class <- allele.diag(filetofix)

alleles.nochange  <- flip.class == "nochange"
alleles.rev <- flip.class == "reverse"
alleles.comp <- flip.class == "comp"
alleles.revcomp <- flip.class == "revcomp"
alleles.impossible <- flip.class == "impossible" 

filetofix[alleles.nochange, REF:=V4]
filetofix[alleles.rev, REF:=V5]
filetofix[alleles.comp, REF:=switches[V4]]
filetofix[alleles.revcomp, REF:=switches[V5]]

filetofix <- filetofix[!alleles.impossible,]
# Remove G/C and A/T SNPs if there are complemented or revcomp SNPs
if(length(flip.class[alleles.revcomp]) +  length(flip.class[alleles.comp]) > 0){
  	filetofix[,alleles:=paste(REF,ALT,sep="/")]
	alleles.toremove <- c("A/T","T/A","G/C","C/G")
	filetofix  <- filetofix[!alleles %in% alleles.toremove,]
	filetofix[,alleles:=NULL]
}


filetofix[,c("V4", "V5", "coord") := NULL]
filetofix <- unique(filetofix)
message("Processed ", paste(filebasename, "-fixed.tsv.gz", sep="")," has ", nrow(filetofix), " SNPs. Writing...")
fwrite(filetofix, paste(filebasename, "-fixed.tsv.gz", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
