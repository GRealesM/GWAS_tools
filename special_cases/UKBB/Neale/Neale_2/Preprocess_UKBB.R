# Pre-processing UKBB files
# Once we have extracted only SNPs present in e IMD, cell, and cyto manifests, we pre-process them prior to pipelining to transform the "variant" column, containing Chr, bp, ref and alt alleles, to our format. Then recompress and save.

library(data.table)
setDTthreads(0)

i <- commandArgs(trailingOnly = TRUE)
message("Pre-processing ", i)
f <- fread(i, tmpdir="tmp")
f[, c("CHR", "BP", "REF", "ALT") := tstrsplit(variant, ":", fixed=TRUE)][, variant :=NULL]
setcolorder(f, c("CHR","BP", "REF","ALT","minor_allele","minor_AF","expected_case_minor_AC","low_confidence_variant", "n_complete_samples", "AC", "ytx","beta",  "se", "tstat","pval"))
newname <- sub("\\-ft1.tsv.gz", "-rh.tsv.gz", i)
fwrite(f, newname, sep = "\t", quote = FALSE)
message("Reheading done!")

