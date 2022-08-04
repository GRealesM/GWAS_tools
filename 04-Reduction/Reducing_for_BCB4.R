#########################################
## DATASET FILTERING BY THE MANIFEST  ###
#########################################

# Introduction: This code is meant to pre-process our big files (located at 02-Liftovered/), 
# filtering them by the SNPs in the SNP manifest, # to make them more manageable prior to 
# project them onto the Blood cell basis

##############################################
### LOAD LIBRARIES AND SET REQUIRED FUNCTIONS
##############################################

library(data.table)
library(magrittr)
library(IMDtools)
setDTthreads(0)
load("~/rds/rds-cew54-basis/03-Bases/BCB4/basis_data/BCB4-Astle_1e4_p1e8.RData")

#############################
## DATA INPUT
#############################

# Get the range of files from args
# This time we'll use a different strategy, involving array jobs.
args <- commandArgs(trailingOnly = TRUE) # Full path to the file

name <- strsplit(args, "\\/|-")[[1]] %>% .[length(.) - 1]
newname <- paste0(name, "-ft.tsv")
message("Reducing ", name, " in ", args, " for Blood cell basis 4.")
input <- fread(args, tmpdir = "tmp")

mincold <- c("CHR38", "BP38", "REF", "ALT", "BETA", "SE", "P")
input <- input[, ..mincold]
M <- g.align(input, SNP.manifest)
dups <- M[duplicated(pid), pid]
if(length(dups) > 0) M <- M[!pid %in% dups] # Remove dups

fwrite(M, paste0("~/rds/rds-cew54-basis/03-Bases/BCB4/reduced_datasets/",newname), sep = "\t")
cat("Done!\n")

