# Computing BETA and/or SE from Z
# Version 1.0

# This script will use data table, alele frequencies and case/control sizes to compute BETA and SE when they're missing, but Z is available.
# File Requirements
# - To have a Z column
# - To have a ALT_FREQ column containing Alt allele frequency (either from the original data or computed using a reference panel
# 
# Three arguments must be provided (in this order)
# (1) The file name
# (2) Number of controls
# (3) Number of cases
#
#

library(data.table)

args <- commandArgs(trailingOnly = TRUE) 
if(length(args) != 3) stop("Three arguments (File name, N0 and N1) are required, please check")

N0  <- as.numeric(args[2])
N1 <- as.numeric(args[3])

N <- N0+N1 # Sample size
s <- N1/N # Proportion of cases

filebasename <- strsplit(args[1], split="-|.tsv.gz")[[1]][1]
ds <- fread(args[1])

essential <- c("ALT_FREQ","Z")
if(!all(essential %in% names(ds))) stop("There are missing columns, either Z or ALT_FREQ. Please check!")

ds[,SE:= sqrt( 1 / (2 * N * ALT_FREQ * (1 - ALT_FREQ) * s * (1 - s)))][,BETA:=Z*SE]
fwrite(ds, paste0(filebasename,"-withBETAandSE.tsv.gz"), sep="\t")


