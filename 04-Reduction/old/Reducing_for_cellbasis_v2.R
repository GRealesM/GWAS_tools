#########################################
## DATASET FILTERING BY THE MANIFEST
#########################################

# Introduction: This code is meant to pre-process our big files (located at 02-Liftovered/), filtering them by the SNPs in the SNP manifest, to make them more manageable prior to project them onto the Blood cell basis

##############################################
### LOAD LIBRARIES AND SET REQUIRED FUNCTIONS
##############################################

library(data.table)
library(magrittr)
setDTthreads(0)

load("~/rds/rds-cew54-basis/03-Bases/cell_basis_v2/cell-basis-sparse-2.0.RData")
SNP.manifest[,alleles:=paste(ref_a1,ref_a2, sep="/")]

# Load some helper functions
### g.complement, as found in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)
g.complement <- function (x) {
  x <- toupper(x)
  switches <- c(A = "t", T = "a", C = "g", G = "c")
  for (i in seq_along(switches)) x <- sub(names(switches)[i], switches[i], x)
  toupper(x)
}

### g.rev, as found in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)
g.rev <- function (x, sep = "/") {
  sapply(strsplit(x, sep), function(g) paste(rev(g), collapse = "/"))
}

### g.class, as found in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)
g.class <- function (x, y) {
  if (!identical(names(x), names(y))) 
    stop("x and y must relate to same SNPs")
  mat <- matrix(FALSE, length(x), 4, dimnames = list(names(x), c("nochange", "rev", "comp", "revcomp")))
  mat[, "nochange"] <- x == y
  mat[, "rev"] <- x == g.rev(y)
  mat[, "comp"] <- x == g.complement(y)
  mat[, "revcomp"] <- x == g.rev(g.complement(y))
  indels <- x %in% c("I/D", "D/I")
  if (any(indels)) 
    mat[indels, c("comp", "revcomp")] <- FALSE
  ret <- character(nrow(mat))
  rs <- rowSums(mat)
  if (length(wh <- which(rs > 1))) 
    ret[wh] <- "ambig"
  if (length(wh <- which(rs == 0))) 
    ret[wh] <- "impossible"
  if (length(wh <- which(rs == 1))) 
    ret[wh] <- colnames(mat)[apply(mat[wh, , drop = FALSE], 1, which)]
  return(ret)
}


#############################
## DATA INPUT
#############################

# Get the range of files from args
# This time we'll use a different strategy, involving array jobs.
# Since we want to control the size of the batch to be similar in all cases, we'll create another variable to store the index + number of extra jobs
args <- commandArgs(trailingOnly = TRUE)
#if(length(args) > 0) {
#	args <- as.numeric(args)
#	start_idx  <- args
#	end_idx <- start_idx + 99 # Set the amount of files we want to process each time (and set array job acconrdingly). Here we'll process, for example, files 1-99, so next job can start from 101.
#	if(end_idx > length(dir(pattern="*-hg38.tsv.gz"))){
#		end_idx <- length(dir(pattern="*-hg38.tsv.gz")) 
#	}
#	files <- dir(pattern="*-hg38.tsv.gz")[start_idx:end_idx]
#} else{
#	files <- dir(pattern="*-hg38.tsv.gz") 
#}

# We changed the strategy here. If we have many files, we can submit an array job and process each of them individually, we just need to supply the name of the file we want to reduce
if(length(args) == 1){
	message("You supplied a file. Processing that file")
	files  <- args
} else {
	message("Processing all (*-hg38.tsv.gz) files")
        files  <- dir(pattern="*-hg38.tsv.gz") # If we have only a few files, by not supplying any argument we can process them all in one go.
}

for(i in files){

cat("Working on ", i, ".\n", sep = "")
input <- fread(i)
input <- input[, c("SNPID", "CHR38", "BP38","REF","ALT", "BETA", "SE", "P")]
input[,alleles:=paste(REF,ALT,sep="/")][,pid:=paste(CHR38,BP38,sep=":")]

M <- merge(input,SNP.manifest[,.(pid,alleles)], by='pid', suffixes=c("",".manifest"))
M$alleles <- toupper(M$alleles)
M$REF <- toupper(M$REF)
M$ALT <- toupper(M$ALT)

rm(input)

  # Check if everything is alright  
  if(!all(g.class(M$alleles.manifest, M$alleles)== "nochange")){
    allele_diagnostics <- g.class(M$alleles.manifest, M$alleles)
    alleles_to_flip <-  allele_diagnostics == "rev"
    alleles_to_comp <- allele_diagnostics == "comp"
    alleles_to_revcomp <- allele_diagnostics == "revcomp"
    cat("Some SNPs have to be flipped. ", sum(alleles_to_flip), " to flip, ", sum(alleles_to_comp), " to find their complement, and ", sum(alleles_to_revcomp), " to find their reverse complement.\n")
    M$alleles[alleles_to_flip] <- unlist(g.rev(M$alleles[alleles_to_flip]))
    M$alleles[alleles_to_comp] <- g.complement(M$alleles[alleles_to_comp])
    M$alleles[alleles_to_revcomp] <- unlist(g.rev(g.complement(M$alleles[alleles_to_revcomp])))
    M$REF <- sapply(strsplit(M$alleles, "/"), `[`, 1)
    M$ALT <- sapply(strsplit(M$alleles, "/"), `[`, 2)
    M$BETA[alleles_to_flip] <- M$BETA[alleles_to_flip]*-1
    M$BETA[alleles_to_revcomp] <- M$BETA[alleles_to_revcomp]*-1

   # NOTE: I introduced the following bit from milcytokine_basis on to guarantee that we have no ambiguity nor duplicated SNPs
    if(!all(g.class(M$alleles.manifest, M$alleles)== "nochange")){
	M  <- M[g.class(M$alleles.manifest, M$alleles)== "nochange",]
	}

    rm(alleles_to_flip, alleles_to_comp, alleles_to_revcomp)
  }
M[, c("alleles", "alleles.manifest"):=NULL]
M <- unique(M)
newname <- strsplit(i, split = "-")[[1]][1]

fwrite(M, paste0("~/rds/rds-cew54-basis/03-Bases/cell_basis_v2/reduced_datasets/",newname,"-ft.tsv"), sep = "\t")
cat("Done!\n")

}
 
