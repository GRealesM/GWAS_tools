## Undoing scales
## Version 1.0

## We all make mistakes, let's undo "corrected" OR scales



library(data.table)

for(i in dir(pattern="*tsv.gz")){
	message("Undoing scale correction in ", i)
	file  <- fread(i)
	file[,c("BETA","SE"):=list(orig.BETA,orig.SE)]
	file[,c("orig.BETA","orig.SE"):=NULL]
	fwrite(file, i, sep="\t")
}
