#########################################
## DATASET FILTERING BY THE MANIFEST
#########################################

# Introduction: This code is meant to pre-process our big files (located at 04-Liftovered/), filtering them by the SNPs in the SNP manifest, to make them more manageable prior to project them onto the IMDbasis

##############################################
### LOAD LIBRARIES AND SET REQUIRED FUNCTIONS
##############################################


library(cupcake)
library(data.table)
library(magrittr)


# Since SNP.manifest at Cupcake is in hg19, we liftover it to hg38
if(!"SNP.manifest.38.tsv" %in% dir()){
SNP.manifest <- copy(cupcake::SNP.manifest)
manifest.build.tanslator <- read.table("Manifest_build_translator.tsv", header = TRUE, sep = "\t")
SNP.manifest <- data.table(pid38 = paste(manifest.build.tanslator$CHR38, manifest.build.tanslator$BP38, sep = ":"), SNP.manifest)
SNP.manifest[,pid:=NULL]
SNP.manifest[,alleles:=paste(ref_a1,ref_a2,sep="/")]
write.table(SNP.manifest, "SNP.manifest.38.tsv", quote = F, sep = "\t", row.names = F)
}

SNP.manifest <- fread("SNP.manifest.38.tsv")

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

for(i in dir(pattern = "*gz")){

#for(i in dir("../04-Liftovered/", pattern = "*gz")[21:length(dir("../04-Liftovered/", pattern = "*gz"))]){
cat("Working on ", i, ".\n", sep = "")
input <- fread(i)
input <- input[, c("SNPID", "CHR38", "BP38","REF","ALT", "BETA", "SE", "P")]
input[,alleles:=paste(REF,ALT,sep="/")][,pid38:=paste(CHR38,BP38,sep=":")]

M <- merge(input,SNP.manifest[,.(pid38,alleles)], by='pid38', suffixes=c("",".manifest"))
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
newname <- strsplit(i, split = "-")[[1]][1]

write.table(M, paste("../IMD_basis/reduced_datasets/",newname,"-ft.tsv", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
cat("Done!\n")

}
 
