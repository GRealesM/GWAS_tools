# Computing BETA and/or SE from Z (imputted data version).
# Version 1.0

# This script will use data table, alele frequencies, effective size, and case/control sizes to compute BETA and SE from imputted datasets (eg. from ssimp) when they're missing, but Z is available.
# File Requirements
# - A Z column
# - An ALT_FREQ column containing Alt allele frequency (either from the original data or computed using a reference panel
# - An N .imp(effective sample size) column, usually provided by the imputation program
# 
# Three arguments must be provided (in this order)
# (1) The file name
# (2) Number of controls
# (3) Number of cases
# 
library(data.table)
# Get all help that you can
setDTthreads(0)

args <- commandArgs(trailingOnly = TRUE) 
if(length(args) != 3) stop("Three arguments (File name, N0 and N1) are required, please check")

N0  <- as.numeric(args[2])
N1 <- as.numeric(args[3])

N <- N0+N1 # Sample size
s <- N1/N # Proportion of cases

filebasename <- strsplit(args[1], split="-|.tsv.gz")[[1]][1]


f.z2b <- function(z, af, n, s) {
    # z = imputed z statistics
    # af = allele frequency
    # n = sample size (effective)
    # s = fraction of samples that are cases
    se.b <- 1/sqrt(2* af * (1-af) * n * s * (1-s))
    b <- z * se.b
    return(list(b, se.b))
}

ds <- fread(args[1])

essential <- c("ALT_FREQ","Z", "N.imp")
if(!all(essential %in% names(ds))) stop("There are missing columns, either Z or ALT_FREQ. Please check!")

ds  <- ds[N.imp > 0] # Remove SNPs with 0 effective sample sizes
ds[,c("BETA","SE"):=f.z2b(z=Z, af=ALT_FREQ, n=N.imp,s=s)]

fwrite(ds, paste0(filebasename, "-withBETAandSE.tsv.gz"), sep="\t", quote=FALSE, na="NA")




