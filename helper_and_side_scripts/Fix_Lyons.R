# Fixing Lyons EGPA
# Goal: To try and adjust BETA in Lyons EGPA files

library(data.table)

convertORscale <- function(x,cp) x/(cp * (1-cp))


EGPA <- fread("../97-Archive/EGPA_Lyons_31719529_1.tsv.gz")
EGPAneg <- fread("../97-Archive/EGPAANCAn_Lyons_31719529_1.tsv.gz")
EGPAMPO <- fread("../97-Archive/EGPAMPO_Lyons_31719529_1.tsv.gz")

cpEGPA <- 534 / (6688 + 535)
cpEGPAneg <- 352 / (6688 + 352)
cpEGPAMPO <- 159 / (6688 + 159)

EGPA$convertedBETA <- convertORscale(EGPA$BETA, cpEGPA)
EGPAneg$convertedBETA <- convertORscale(EGPAneg$BETA, cpEGPAneg)
EGPAMPO$convertedBETA <- convertORscale(EGPAMPO$BETA, cpEGPAMPO)
EGPA$convertedSE <- convertORscale(EGPA$SE, cpEGPA)
EGPAneg$convertedSE <- convertORscale(EGPAneg$SE, cpEGPAneg)
EGPAMPO$convertedSE <- convertORscale(EGPAMPO$SE, cpEGPAMPO)
EGPA$convertedOR <- exp(EGPA$convertedBETA)
EGPAneg$convertedOR <- exp(EGPAneg$convertedBETA)
EGPAMPO$convertedOR <- exp(EGPAMPO$convertedBETA)

EGPA[, `:=`(BETA=NULL, SE=NULL, id=NULL)]
EGPAneg[, `:=`(BETA=NULL, SE=NULL, id=NULL)]
EGPAMPO[, `:=`(BETA=NULL, SE=NULL, id=NULL)]
setnames(EGPA, old = c("convertedBETA", "convertedSE", "convertedOR"), new=c("BETA", "SE", "OR"))
setnames(EGPAneg, old = c("convertedBETA", "convertedSE", "convertedOR"), new=c("BETA", "SE", "OR"))
setnames(EGPAMPO, old = c("convertedBETA", "convertedSE", "convertedOR"), new=c("BETA", "SE", "OR"))

fwrite(EGPA, "EGPA_Lyons_31719529_1.tsv.gz", row.names=FALSE, quote=FALSE, sep="\t")
fwrite(EGPAneg, "EGPAANCAn_Lyons_31719529_1.tsv.gz", row.names=FALSE, quote=FALSE, sep="\t")
fwrite(EGPAMPO, "EGPAMPO_Lyons_31719529_1.tsv.gz", row.names=FALSE, quote=FALSE, sep="\t")
