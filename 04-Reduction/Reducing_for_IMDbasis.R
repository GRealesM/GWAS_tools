#########################################
## DATASET FILTERING BY THE IMD MANIFEST
#########################################

# Introduction: This code is meant to pre-process our big files (located at 02-Liftovered/), filtering them by the SNPs in the SNP manifest, to make them more manageable prior to project them onto the IMDbasis

##############################################
### LOAD LIBRARIES AND SET REQUIRED FUNCTIONS
##############################################

library(cupcake)
library(data.table)
library(magrittr)
library(IMDtools)
setDTthreads(0)


# Since SNP.manifest at Cupcake is in hg19, we liftover it to hg38
if(!"SNP.manifest.38.tsv" %in% dir("~/rds/rds-cew54-basis/GWAS_tools/04-Reduction/")){
SNP.manifest <- copy(cupcake::SNP.manifest)
manifest.build.tanslator <- read.table("~/rds/rds-cew54-basis/GWAS_tools/01-Pipeline/Manifest_build_translator.tsv", header = TRUE, sep = "\t")
SNP.manifest <- data.table(pid38 = paste(manifest.build.tanslator$CHR38, manifest.build.tanslator$BP38, sep = ":"), SNP.manifest)
SNP.manifest[,pid:=NULL]
SNP.manifest[,alleles:=paste(ref_a1,ref_a2,sep="/")]
write.table(SNP.manifest, "~/rds/rds-cew54-basis/GWAS_tools/04-Reduction/SNP.manifest.38.tsv", quote = F, sep = "\t", row.names = F)
}
SNP.manifest <- fread("~/rds/rds-cew54-basis/GWAS_tools/04-Reduction/SNP.manifest.38.tsv")

#######################
## Prepare manifest  ##
#######################

SNP.manifest[, c("CHR38", "BP38") := tstrsplit(pid38, split = ":")]
setnames(SNP.manifest, c("ref_a1", "ref_a2"), c("REF", "ALT"))
SNP.manifest <- SNP.manifest[, .(CHR38, BP38, REF, ALT)]


#############################
## DATA INPUT
#############################
# Get the range of files from args
# This time we'll use a different strategy, involving array jobs.
args <- commandArgs(trailingOnly = TRUE) # Full path to the file

# We changed the strategy here. If we have many files, we can submit an array job and process each of them individually, we just need to supply the name of the file we want to reduce
if(length(args) == 1){
	message("You supplied a file. Processing that file: ", args)
	files  <- args
} else {
	message("Processing all (*-hg38.tsv.gz) files")
        files  <- dir(pattern="*-hg38.tsv.gz") # If we have only a few files, by not supplying any argument we can process them all in one go.
}


for(i in files){

cpath <- strsplit(i, split ="\\/")[[1]] # In case a path is supplied
cpath <- cpath[length(cpath)]

name <- strsplit(cpath, split = "-")[[1]][1]
newname <- paste0(name, "-ft.tsv")
message("Reducing ", name, " for IMD basis.")
input <- fread(i, tmpdir = "tmp")

mincold <- c("CHR38", "BP38", "REF", "ALT", "BETA", "SE", "P")
input <- input[, ..mincold]
M <- g.align(input, SNP.manifest)
dups <- M[duplicated(pid), pid]
if(length(dups) > 0) M <- M[!pid %in% dups] # Remove dups
setnames(M, "pid", "pid38") # We need pid column to be named pid38

fwrite(M, paste0("~/rds/rds-cew54-basis/03-Bases/IMD_basis/reduced_datasets/",newname), sep = "\t")
cat("Done!\n")

}
 