# Script for merging genomic positions and meta-analysis files

args <- commandArgs(trailingOnly = TRUE)

library(data.table)

genpos <- fread(args[1], na.strings="")
meta <- fread(args[2])

traitname <- strsplit(args[1], split = "_")[[1]][1]

fullmeta <- merge(genpos, meta, by="SNPID", all.y=TRUE)
fullmeta[,REF:=toupper(REF)][,ALT:=toupper(ALT)]

fwrite(fullmeta, paste(traitname, "_Reales_up_1-hg38.tsv.gz", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t", na = "NA")
