#########################################
## DATASET FILTERING BY THE MANIFEST
#########################################

# Introduction: This code is intended to reduce all files in cytokine_basis filtering them by the SNPs in the cytokine13_basis SNP manifest, since cytokin13 is a subset of cytokine_basis.

##############################################
### LOAD LIBRARIES AND SET REQUIRED FUNCTIONS
##############################################

library(data.table)
library(magrittr)

# # Get the range of files from args
# args <- commandArgs(trailingOnly = TRUE)
# args <- as.numeric(args)

load("~/rds/rds-cew54-basis/GWAS_tools/03-Bases/cytokine13_basis/cytokine13-basis-sparse-1.0.RData")
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

for(i in dir("~/rds/rds-cew54-basis/03-Bases/cytokine_basis/reduced_datasets/", pattern="-ft.tsv")){

cat("Working on ", i, ".\n", sep = "")
filepath  <- file.path("../03-Bases/cytokine_basis/reduced_datasets/", i)
input <- fread(filepath)
input <- input[, c("SNPID", "CHR38", "BP38","REF","ALT", "BETA", "SE", "P")]
input[,alleles:=paste(REF,ALT,sep="/")][,pid:=paste(CHR38,BP38,sep=":")]

# I included the next line from milmetacytokine on, with the hope that it would make processing faster.
input  <- input[pid %in% SNP.manifest$pid,]

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

write.table(M, paste("~/rds/rds-cew54-basis/GWAS_tools/03-Bases/cytokine13_basis/reduced_datasets/",newname,"-ft.tsv", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
cat("Done!\n")

}
 
