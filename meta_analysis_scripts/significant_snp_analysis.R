#################################
### SIGNIFICANT SNP ANALYSIS ####
#################################

# Background: This script is intended to analyse the significant SNPs found in the Meta-analysis and its source files, with the following goals:
# (1) Are there new hits after meta-analysis?
# (2) For significant SNPs in source files, did they become more significant in the Meta-analysis, or not?

# Loading libraries
library(bigreadr)
library(data.table)


# Initial setting
thres  <- 5e-8
new_df  <- data.table()
old_df  <- data.table()
lost_df <- data.table()

# We take advantage that files in meta-analysis and input have the hg38 label
files <- dir(path="..", pattern = "hg38.tsv.gz")
traits <- unique(sapply(strsplit(files, "_"), `[`, 1))

for(i in traits){

message("Working on ", i)
files_for_i  <- files[grepl(i,files)]

# In alphabetic order, meta-analysis (Reales) will always be third
metalabel <- files_for_i[3]
input1label <- files_for_i[1]
input2label <- files_for_i[2]

message("Loading corresponding files: ", input1label,", ", input2label,", and ",metalabel,".")  
meta  <- bigreadr::fread2(paste("../",metalabel, sep=""), data.table = TRUE, nThread = 8)
input1 <- bigreadr::fread2(paste("../",input1label, sep=""), data.table = TRUE, nThread = 8)
input2 <- bigreadr::fread2(paste("../",input2label, sep=""), data.table = TRUE, nThread = 8)

# Filter unnecessary columns in inputs
input1  <- input1[, c("SNPID", "CHR38", "BP38", "REF", "ALT", "BETA", "SE", "P")]
input2  <- input2[, c("SNPID", "CHR38", "BP38", "REF", "ALT", "BETA", "SE", "P")]

# Checking for significant SNPs in each file
sigmeta  <- meta[P < thres, SNPID]
if(length(sigmeta) == 0){ message("No significant hits in meta-analysis")
	} else { message(length(sigmeta), " significant hits in meta-analysis")
}
siginput1  <- input1[P < thres, SNPID]
if(length(siginput1) == 0){ message("No significant hits in input1")
	} else { message(length(siginput1), " significant hits in input1")}
siginput2  <- input2[P < thres, SNPID]
if(length(siginput2) == 0){ message("No significant hits in input2")
	} else { message(length(siginput2), " significant hits in input2")}
sigtotal <- unique(c(sigmeta,siginput1,siginput2)) # Significant SNPs in at least one dataset
siginputtotal  <- unique(c(siginput1,siginput2)) # Significant SNPs in at least one of the inputs
newhits <- setdiff(sigmeta,siginputtotal) # Significant SNPs only in meta
if(length(newhits) == 0){ 
	message("No new significant hits for ",i," in meta-analysis! :(")
} else { message(length(newhits), " new significant hits in ", i," in meta-analysis! Great!")}

# Subset data.tables to make tables

# First, for new hits
if(length(newhits > 0)){
	newhits_df  <- rbindlist(list(meta[SNPID %in% newhits, file:=metalabel],
				      input1[SNPID %in% newhits, file:=input1label],
				      input2[SNPID %in% newhits, file:=input2label]), fill = TRUE)
	newhits_df[,trait:=i]
	new_df  <- rbind(new_df,newhits_df, fill=TRUE)
} else { message("Nothing to add to new_df at this point.")}

# Second, to compare how old hits changed in the meta-analysis
oldhits <- intersect(sigmeta, siginputtotal) # Significant SNPs in meta and in at least one of the inputs
if(length(oldhits > 0)){
	message("There are ", length(oldhits) ," significant hits in meta that were significant in at least one of the inputs")
	oldhits_df  <- rbindlist(list(meta[SNPID %in% sigmeta, file:=metalabel],
				      input1[SNPID %in% sigmeta, file:=input1label],
				      input2[SNPID %in% sigmeta, file:=input2label]), fill=TRUE)
	oldhits_df[,trait:=i]
	old_df  <- rbind(old_df,oldhits_df, fill=TRUE)
} else { message("Nothing to add to old_df at this point.")}

# Third, to check how lost hits (sig in input, not sig in meta) changed
losthits  <- setdiff(siginputtotal, sigmeta)
if(length(losthits > 0)){
	losthits_df  <- rbindlist(list(meta[SNPID %in% losthits, file:=metalabel],
				      input1[SNPID %in% losthits, file:=input1label],
				      input2[SNPID %in% losthits, file:=input2label]), fill=TRUE)
	losthits_df[,trait:=i]
	lost_df  <- rbind(lost_df, losthits_df, fill=TRUE)
} else { message("Nothing to add to lost_df at this point.")}

rm(meta,input2,input2)

}

bigreadr::fwrite2(new_df, "new_hits.tsv", sep="\t", row.names=F, nThread=1)
bigreadr::fwrite2(old_df, "old_hits.tsv", sep="\t", row.names=F, nThread=1)
bigreadr::fwrite2(lost_df, "lost_hits.tsv", sep="\t", row.names=F, nThread=1)







