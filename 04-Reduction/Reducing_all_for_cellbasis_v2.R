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

# Get the range of files from args
# This time we'll use a different strategy, involving array jobs.
# Since we want to control the size of the batch to be similar in all cases, we'll create another variable to store the index + number of extra jobs
args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)


start_idx  <- args
end_idx <- start_idx + 9
if(end_idx > length(dir("../../02-Processed/", pattern="*tsv.gz"))){
	end_idx <- length(dir("../02-Processed/", pattern="*tsv.gz")) 
}

load("../../03-Bases/cell_basis_v2/cell-basis-sparse-2.0.RData")
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

for(i in dir("../../02-Processed/", pattern="*tsv.gz")[start_idx:end_idx]){

cat("Working on ", i, ".\n", sep = "")
filepath  <- file.path("../../02-Processed/", i)
input <- fread(filepath)
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
    
    rm(alleles_to_flip, alleles_to_comp, alleles_to_revcomp)
  }
M[, c("alleles", "alleles.manifest"):=NULL]
M <- unique(M)

# Sometimes more than one SNP (with different alleles) are mapped to the same position. To prevent SNPs in the same position, but not in the manifest to proceed, let's do a quick check.
# This way only SNPs with exactly the same pid, and both alleles (and aligned) as the SNP.manifest will proceed.

if(any(duplicated(M$pid))){
 snp.M <- paste(M$pid, M$REF, M$ALT, sep = ":")
 snp.man <- paste(SNP.manifest$pid, SNP.manifest$ref_a1, SNP.manifest$ref_a2, sep = ":")
 M <- M[snp.M %in% snp.man,]
}

newname <- strsplit(i, split = "-")[[1]][1]

write.table(M, paste("../../03-Bases/cell_basis_v2/reduced_datasets/",newname,"-ft.tsv", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
cat("Done!\n")

}
 
