######################################
### PRE-PROCESSING SAKAUE FILES  #####
######################################

# This script is meant to address the 215 datasets from Sakaue et al. 2021, corresponding to GWAS on 215 traits in Japanese subjects from Biobank Japan.
# Many files come with autosomes and X chromosome in separate files, so we'll combine them into one file.
# The traits include diseases, biomarkers, and drug administration. Some of these traits used BOLT-LMM, so we need to convert from linear to OR scale.
# Also, there's a manifest of the columns, so we'll use it to prepare the columns for pipeline treatment.

# Use arguments supplied by array job
args <- as.numeric(commandArgs(trailingOnly=TRUE))

# Load libraries
library(data.table)
library(IMDtools)
setDTthreads(8)

dirs <- dir(pattern="hum0197*") # All directories
indir <- dir(path=dirs[args]) # Files in directory
namedict <- fread("namedict")

# Merge if two files are provided.

if(length(indir) == 2){ # In case there are two files
	message("Merging...")
	auto <- fread(paste(dirs[args], indir[1], sep="/")) 
	x <-  fread(paste(dirs[args], indir[2], sep="/"))
	ds <- rbindlist(list(auto,x), use.names=TRUE, fill=TRUE)
}else{
	ds <- fread(paste(dirs[args], indir, sep="/"))

}

# Detect if it's a biomarker, and apply sdY correction

if("P_BOLT_LMM_INF" %in% names(ds)){ # Check if quantitative trait, treated with BOLT-LMM
	message("Adjusting sdY...")
	n <- namedict[original.name == dirs[args]]$n
	ds[,c("original.BETA","original.SE"):=list(BETA,SE)][, c("BETA","SE"):=sdY.correction(BETA,SE,A1FREQ,n)]

}

# Lastly, rename
setnames(ds, c("ALLELE1","ALLELE0","A1FREQ","P_BOLT_LMM_INF", "Allele1","Allele2","AF_Allele2","p.value"),
	     c("ALT","REF","ALT_FREQ","P","REF","ALT", "ALT_FREQ","P"), skip_absent=TRUE)

# Save the file
newname <- namedict[original.name == dirs[args]]$new.name
fwrite(ds, newname, sep="\t")
