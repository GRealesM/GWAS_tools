# Fetching coordinates
# version 1

# This script will use 1000 Genomes Phase III to retrieve genomic coordinates from SNPID

library(data.table)
setDTthreads(10)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 1) stop("Only one file per run at the time. Sorry!")

message("Loading ",args[1]," dataset.")

ds <- fread(args[1])
if(!"SNPID" %in% names(ds)) stop("SNPID column is required! Please check!")

reference.path <- "~/rds/rds-cew54-basis/95-1000genomes/reference_hg19/"
filebasename <- strsplit(args[1], split="-|.tsv.gz")[[1]][1]
snplist <- ds$SNPID


message("Extracting coordinates from reference. This may take a while.")
coordlist <- lapply(c(1:22, "X", "Y"), function(i){
	message("Working on chr ",i,".")
	bim <- fread(paste0(reference.path,"chr",i,".bim"), select=c("V1","V4","V2"),col.names=c("CHR","BP", "SNPID"))
	bim <- bim[SNPID != ".",]
	multrsids  <- bim[grep(";", bim$SNPID),]
	
	if(nrow(multrsids) > 0 ){
	bim  <- bim[!grep(";", bim$SNPID),]
	multrsids <- multrsids[,.(CHR, BP, SNPID, splitted=unlist(strsplit(SNPID, ";"))), by=seq_len(nrow(multrsids))][,SNPID:=splitted][,c("seq_len", "splitted"):=NULL]
	bim  <- rbind(bim,multrsids)
	}
	idx <- bim$SNPID %in% snplist
	bim <- bim[idx,]
	message(nrow(bim)," variants were found in chr ", i,"!")
	bim
})
coords <- rbindlist(coordlist)
out <- merge(ds, coords, by="SNPID")
out <- out[order(CHR,BP),]

fwrite(out, paste0(filebasename, "-withcoords.tsv.gz"), sep="\t")


