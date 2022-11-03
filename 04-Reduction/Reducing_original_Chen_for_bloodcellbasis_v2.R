#############################################################
## DATASET FILTERING BY THE BLOOD CELL BASIS v2 MANIFEST  ###
#############################################################

# Some Chen datasets had weird projections. We performed a PCA on reduced datasets (sdY-transformed) BETAs to see how raw data behaved. 
# This confirmed that East Asian datasets didn't group with the rest, for whatever reason.
# Previously, by investigating HBC_6, with weird projection values (see `~/rds/rds-cew54-basis/GWAS_tools/special_cases/exploring_HBC6.R`) 
# we realised that this problem was probably caused by some SNPs that had too small SE, which resulted in absurd values when sdY-correction was applied.
# Here we reduced and aligned Chen datasets to the manifest, using original BETA (labelled "orig.BETA") to see how much of these issues correspond to the 
# sdY transformation.

##############################################
### LOAD LIBRARIES AND SET REQUIRED FUNCTIONS
##############################################

library(data.table)
library(magrittr)
library(IMDtools)
setDTthreads(0)
load("~/rds/rds-cew54-basis/03-Bases/cell_basis_v2/cell-basis-sparse-2.0.RData")

#######################
## Prepare manifest  ##
#######################

SNP.manifest[, c("CHR38", "BP38") := tstrsplit(pid, split = ":")]
setnames(SNP.manifest, c("ref_a1", "ref_a2"), c("REF", "ALT"))
SNP.manifest <- SNP.manifest[, .(CHR38, BP38, REF, ALT)]


#############################
## DATA INPUT
#############################

# Note the procedure here is a bit different to the other "Reducing" scripts
# We'll provide a number (1 to 89) corresponding to an index for the 89 Chen datasets.
# Then we'll define internally the list of Chen datasets
args <- commandArgs(trailingOnly = TRUE) # Number (1 to 89)
args <- as.numeric(args)

ptorawfiles <- "~/rds/rds-cew54-basis/02-Processed/"
chenrawfiles <- grep("Chen_", dir(ptorawfiles), value = TRUE)

name <- chenrawfiles[args] # Select file
basename  <- gsub("-hg38.tsv.gz", "", name)    # Remove hg38.. extension
newname <- paste0(basename, "-ft.tsv")         # Output name

message("Reducing original Chen dataset ", basename, " for Blood cell basis 2.")
input <- fread(paste0(ptorawfiles, name), tmpdir = "tmp")

mincold <- c("CHR38", "BP38", "REF", "ALT", "orig.BETA", "orig.SE", "P") # Important step! Keep original values
input <- input[, ..mincold]
setnames(input, c("orig.BETA", "orig.SE"), c("BETA", "SE"))
M <- g.align(input, SNP.manifest)
dups <- M[duplicated(pid), pid]
if(length(dups) > 0) M <- M[!pid %in% dups] # Remove dups

fwrite(M, paste0("~/rds/rds-cew54-basis/03-Bases/cell_basis_v2/reduced_orig_chen/",newname), sep = "\t") # New directory
cat("Done!\n")

