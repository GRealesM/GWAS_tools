#############################################################
# Scaling BETA and SE from linear scale (BOLT-LMM) or by sdY
#############################################################

# This script will convert BETA and SE from linear to log-odds scale in case-control datasets generated in BOLT-LMM. Alternatively, if the trait is quantitative (ie., N1 = 0), it will estimate sdY and adjust BETA and SE by it.
# It will by default process all ".tsv.gz" files in the directory, and output compressed "-bscorr.tsv.gz" files with old Beta and SE columns named "orig.BETA" and "orig.SE", and new columns named "BETA" and "SE", respectively.

library(data.table)
library(IMDtools)
# devtools::install_github("tidyverse/googlesheets4")
library(googlesheets4)

# Get all help that you can
setDTthreads(0)

# Since we need extra information (N0 and N1), we'll obtain it from our Main_table
gs4_deauth()
message("Fetching metadata...")
metadata <- read_sheet("https://docs.google.com/spreadsheets/d/16B4ANehcS4psdAFReTBQUJLYXuf5_RrpmjSPaASa2Nw/edit?usp=sharing", sheet = 1)[,c("File_ID", "N0", "N1")]
metadata <- data.table(metadata)

for (i in dir(pattern= "*.tsv.gz")){
file <- fread(i)
filemetadata <- metadata[File_ID == i,]
if(nrow(filemetadata) != 1) stop("File couldn't be matched to metadata. Please check.")
N0 <- filemetadata$N0
N1 <- filemetadata$N1
if(!is.na(N1) && N1 > 0){
  message(i," seems to be a case-contol study, scaling...")
  file[,orig.BETA:=BETA][,orig.SE:=SE][,c("BETA", "SE"):=linORscale(beta = BETA, se = SE,N0 = N0, N1 = N1)]
} else {
  message(i, " seems to be a quantitative trait, computing sdY and transforming...")
  file[,orig.BETA:=BETA][,orig.SE:=SE][,c("BETA", "SE"):=sdY.correction(beta = BETA, se = SE, maf = ALT_FREQ, n = N0)]
}
newname <- paste(strsplit(i, split = "-")[[1]][1], "-bsadj.tsv.gz", sep="")
fwrite(file, newname, sep = "\t", quote = FALSE, row.names = FALSE, na="NA")
message("Done!")

}
