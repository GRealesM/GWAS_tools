# Merge freqs back with its original file

#Load libraries
library(data.table)


# Get the range of files from args
args <- commandArgs(trailingOnly = TRUE)

# Helper function


switches <- c(A = "T", T = "A", C = "G", G = "C")

allele.diag.2alleles <- function(x){

	flip.class <- ifelse(x$ALT == x$ALT.1KG & x$REF == x$REF.1KG, 
			     "nochange",
			     ifelse(x$ALT == x$REF.1KG & x$REF == x$ALT.1KG, 
				    "reverse",
				    ifelse(switches[x$ALT]== x$ALT.1KG & switches[x$REF] == x$REF.1KG, 
					   "comp", 
					   ifelse(switches[x$ALT] == x$REF.1KG & switches[x$REF] == x$ALT.1KG,
						  "revcomp",
						  "impossible"))))
	return(flip.class)
}

# Load files

file <- fread(args[1])
freqs <- fread(args[2])

# Prepare and merge

file[, pid:=paste(CHR19,BP19, sep=":")][,c("REF","ALT"):=list(toupper(REF),toupper(ALT))]
freqs[, pid:=paste(CHR19,BP19, sep=":")]

merged  <- merge.data.table(file, freqs[,c("pid", "REF", "ALT", "ALT_FREQ")], by = "pid", all.x=TRUE, suffix = c("", ".1KG"))
message(length(merged$ALT_FREQ[is.na(merged$ALT_FREQ)]), " SNPs couldn't have their allele frequencies computed.") 
# Remove SNPs without freqs
merged  <- na.omit(merged, cols = "ALT_FREQ")

# Let's check if all alleles match. Unlikely unless the dataset is very small and ideal
if(!all(merged$ALT == merged$ALT.1KG)){
	# Sometimes there are trialellic and tetra-allelic SNPs. Those are hard to make a decision on, so we'll remove them
	merged <- merged[!duplicated(pid),]
	# Let's remove all fixed alleles
	merged <- merged[!(ALT_FREQ == 0 | ALT_FREQ == 1), ]

	flip.class <- allele.diag.2alleles(merged)
	
	alleles.rev <- flip.class == "reverse"
	alleles.comp <- flip.class == "comp"
	alleles.revcomp <- flip.class == "revcomp"
	alleles.impossible <- flip.class == "impossible" 

	# Flip ALT_FREQ when necessary
	merged[alleles.rev, ALT_FREQ:=1-ALT_FREQ]
	merged[alleles.revcomp, ALT_FREQ:=1-ALT_FREQ]
	# Remove impossibles 
	merged <- merged[!alleles.impossible,]
	
	# Remove G/C and A/T SNPs if there are complemented or revcomp SNPs
	if(length(flip.class[alleles.revcomp]) +  length(flip.class[alleles.comp]) > 0){
	  	merged[,alleles:=paste(REF,ALT,sep="/")]
		alleles.toremove <- c("A/T","T/A","G/C","C/G")
		merged  <- merged[!alleles %in% alleles.toremove,]
		merged[,alleles:=NULL]
	}
}

merged[,c("REF.1KG", "ALT.1KG","pid"):=NULL]

basename  <- strsplit(args[1], split=c("-|\\.tsv"))[[1]][1]
fwrite(merged, paste(basename, "-withfreqs.tsv.gz", sep=""), sep="\t")

