####### Calculate if SNP rs2476601 BETA and SE is proportional to N0+N1/N0*N1 in Case-control files

library(data.table)
library(ggplot2)
library(ggrepel)

casecontrols <- fread("CaseControl_N0N1_forQC_20200220.tsv")
casecontrols$N0N1Calc <- (casecontrols$N0 + casecontrols$N1)/(as.numeric(casecontrols$N0) * as.numeric(casecontrols$N1))
BETA <- rep(NA, nrow(casecontrols))
SE <- rep(NA, nrow(casecontrols))

for (i in 1:nrow(casecontrols)){
  file <- fread(paste("Filtered_datasets/", casecontrols$File_ID[i], sep =""))
  file <- file[file$pid38 == "1:113834946", ]
  if(nrow(file) == 1) {
    BETA[i] <- file$BETA[1]
    SE[i] <- file$SE[1]
  }
}

casecontrols$BETA <- BETA
casecontrols$SE <- SE
casecontrols$Z <- casecontrols$BETA/casecontrols$SE

seplot <- ggplot(casecontrols, aes(x = log(N0N1Calc), y = log(SE), colour = BETA, label = File_ID)) +
    geom_point()+
    geom_text_repel(data = subset(casecontrols, BETA < -2.5))+
    theme_classic()
seplot

zplot <- ggplot(casecontrols, aes(x = BETA, y = abs(Z), label = File_ID)) +
  geom_point()+
  geom_text_repel(data = subset(casecontrols, Z > 8 | Z < -8))+
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red")+
  theme_classic()
zplot

plot(log(casecontrols$N0N1Calc), casecontrols$BETA)
plot(log(casecontrols$N0N1Calc), log(casecontrols$SE))
