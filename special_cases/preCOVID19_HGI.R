# Pre-treating COVID19 datasets
library(data.table)

datasets <- dir(pattern="-raw.tsv.gz")
setDTthreads(8)
for(i in datasets){
 message("Working on ",i)
 ds <- fread(i)
 setnames(ds, old=c("#CHR", "POS","SNP","all_inv_var_meta_beta", "all_inv_var_meta_sebeta","all_inv_var_meta_p","rsid"), new=c("CHR","BP","SNPalternative", "BETA","SE", "P","SNPID"))
 fwrite(ds, i, sep="\t")

}
