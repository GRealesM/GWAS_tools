# Updating genomic coordinates
# Version 1.0

# Some files come with zero-based genomic positions. This script will correct that.
# This script will use data table with, at least, genomic coordinates.
# File Requirements
# - Have a position column
# 
# Three arguments must be provided (in this order)
# (1) The file name
# (2) Name of position column
# (3) How many bases to move? Use 1 if your dataset is zero-based and you want to add +1 to all positions
#
#

library(data.table)

args <- commandArgs(trailingOnly = TRUE) 
if(length(args) != 3) stop("Three arguments (File name, Position column name and positions to add) are required, please check")

poscol <- args[2]
newpos  <- as.numeric(args[3])


filebasename <- strsplit(args[1], split="-|.tsv.gz")[[1]][1]
ds <- fread(args[1])
if(!poscol %in% names(ds)) stop("Supplied column is not in the header. Please check!")
if(!is.numeric(ds[,get(poscol)])) {
	ds[,get(poscol)] <- as.numeric(ds[,get(poscol)])
}
message("Updating positions!")
ds[,c(poscol):=get(poscol) + newpos]

fwrite(ds, paste0(filebasename,"-updatedpos.tsv.gz"), sep="\t")


