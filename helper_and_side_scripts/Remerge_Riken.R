## REMERGING Japanese ENcyclopedia of GEnetic associations by Riken FILES ####

# This script is intended to merge files from Riken, with different columns

args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

library(data.table)
setDTthreads(threads = 0)



auto <- fread(args[1])
chrx <- fread(args[2])
auto$CHR <- as.character(auto$CHR)
merged <- rbindlist(list(auto,chrx), fill = TRUE)
write.table(merged, args[3], quote = F, row.names = F, sep = "\t")
