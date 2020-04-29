# Selecting significant SNPs 

# Goal: Filter SNPs by their P-values accross all reduced SNPs, and count their occurrence.
# This way we'll get a list of important SNPs for many/most traits, that we can use for QC

library(data.table)

main_df <- data.frame()

for(i in dir("Filtered_datasets/")){
  x <- fread(paste("Filtered_datasets/",i, sep = ""))
  x$filename <- i
  main_df <- rbind(main_df, x, fill= T)
}
filtered_df <- main_df[main_df$P < 1e-5, ]
mainSNPS <- as.data.table(table(filtered_df$SNPID))
setorder(mainSNPS, -N)


# Some BETA/SE/P QC
# I realised that some BETA, SE and P were missing, lets look into that
BETANA <- main_df[is.na(BETA),]
SENA <- main_df[is.na(SE),]
PNA <- main_df[is.na(P),]
