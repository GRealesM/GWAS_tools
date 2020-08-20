# Merge freqs back with its original file

#Load libraries
library(data.table)


# Get the range of files from args
args <- commandArgs(trailingOnly = TRUE)

# Load files

file <- fread(args[1])
freqs <- fread(args[2])

# Prepare and merge

file[, pid:=paste(CHR19,BP19, sep=":")]
freqs[, pid:=paste(CHR19,BP19, sep=":")]

merged  <- merge.data.table(file, freqs[,c("pid", "REF", "ALT", "ALT_FREQ")], by = "pid", all.x=TRUE, suffix = c("", ".1KG"))
message(length(merged$ALT_FREQ[is.na(merged$ALT_FREQ)]), " SNPs couldn't have their allele frequencies computed.") 

# To be conservative, we'll remove all freqs in which the ALT allele in 1000 genomes and in the file don't match
if(!all(merged$ALT == merged$ALT.1KG)){
   message(length(which(merged$ALT != merged$ALT.1KG)), " SNPs had discordant ALT. If this is important, please modify the code to remove/flip them.")
  # merged[merged$ALT != merged$ALT.1KG,  ]$ALT_FREQ  <- NA
}

merged[,c("REF.1KG", "ALT.1KG","pid"):=NULL]

basename  <- strsplit(args[1], split=c("-", ".tsv"))[[1]][1]
fwrite(merged, paste(basename, "-withfreqs.tsv.gz", sep=""), sep="\t")

