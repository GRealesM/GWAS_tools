library(data.table)
setwd("Projecting_dbs/Filtered_datasets/")
CDlange <- fread("CD_DeLange_28067908_1-ft.tsv")
CDbrant <- fread("CD_Brant_27693347_1-ft.tsv")
CDcomb <- merge(CDlange, CDbrant, by = "pid38", suffixes = c(".lange", ".brant"))
plot(CDcomb$BETA.lange, CDcomb$BETA.brant)
abline(lm(BETA.brant ~ BETA.lange, data = CDcomb))

