library(data.table)
setDTthreads(10)



gwama_files <- dir(pattern = "*_[1-5]-raw.tsv.gz")
#specials <- c("PLAV_Chen_32888493_3-raw.tsv.gz", "PLAV_Chen_32888493_4-raw.tsv.gz")
#gwama_files <- setdiff(gwama_files, specials)
mrmegafiles <- dir(pattern = "*_6-raw.tsv.gz")
#
for(i in gwama_files){
filebasename <- strsplit(i, "-")[[1]][1]
message("Working on ",i,".")
ds <- fread(i)
ds[,tmpid:=gsub("_.*","",ds$rs_number)]
ds[,c("CHR","BP"):=tstrsplit(tmpid,":", fixed=TRUE)]
ds <- ds[CHR != "XY",]
ds[,tmpid:=NULL]
setnames(ds, old=c("reference_allele","other_allele","eaf", "beta","se", "z","p-value"), new=c("ALT","REF","ALT_FREQ","BETA","SE","Z","P"))
fwrite(ds, paste0(filebasename,"-hc.tsv.gz"),sep="\t")

}
#
#for(i in specials){
#filebasename <- strsplit(i, "-")[[1]][1]
#message("Working on ",i,".")
#ds <- fread(i)
#ds <- ds[CHR != 24,]
#setnames(ds, old=c("EFFECT_ALLELE","OTHER_ALLELE","EAF","PVAL"), new=c("ALT","REF","ALT_FREQ","P"))
#fwrite(ds, paste0(filebasename,"-hc.tsv.gz"),sep="\t")
#}


for(i in mrmegafiles[6]){
filebasename <- strsplit(i, "-")[[1]][1]
message("Working on ",i,".")
ds <- fread(i)
if("P-value_association" %in% names(ds)) setnames(ds,"P-value_association","P.value_association")
setnames(ds, old=c("MarkerName","Chromosome","Position", "EA","NEA","EAF","beta_0","se_0","P.value_association"), new=c("SNPID","CHR","BP","ALT","REF","ALT_FREQ","BETA","SE","P"))
fwrite(ds, paste0(filebasename,"-hc.tsv.gz"),sep="\t")
}
