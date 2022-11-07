# Process Ishigaki 2022

# The files need some pre-processing before they can be pipelined

library(data.table)
setDTthreads(10)

args  <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)

in.name <- args[1]
out.name <- args[2]

d  <- fread(in.name)
d2 <- d[!grepl("INS", SNP)] # Remove Alu insertions 
d2[, c("CHR", "BP", "REF", "ALT"):=tstrsplit(SNP, "_")]
d2 <- d2[nchar(REF) == 1 & nchar(ALT) == 1] # Keep SNPs only
setnames(d2, c("SNP", "Beta", "Pval"), c("SNPID", "BETA", "P"))
d2 <- d2[, .(SNPID, CHR, BP, REF, ALT, BETA, SE, P)]

fwrite(d2, out.name, sep = "\t")
