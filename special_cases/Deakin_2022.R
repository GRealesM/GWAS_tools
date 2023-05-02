library(data.table)

setwd("/home/gr440/rds/rds-cew54-basis/01-Workshop")
jdm <- fread("JDM_Deakin_35094092_1-raw.tsv.gz")
jdm <- jdm[,.(CHR,BP)]
mf <- fread("../GWAS_tools/01-Pipeline/Manifest_build_translator.tsv")
m18 <- merge(jdm, mf[,.(CHR18, BP18)], by.x = c("CHR", "BP"), by.y= c("CHR18", "BP18"))
m19 <- merge(jdm, mf[,.(CHR19, BP19)], by.x = c("CHR", "BP"), by.y= c("CHR19", "BP19"))
m38 <- merge(jdm, mf[,.(CHR38, BP38)], by.x = c("CHR", "BP"), by.y= c("CHR38", "BP38"))
nrow(m18); nrow(m19); nrow(m38)
# It's hg19!