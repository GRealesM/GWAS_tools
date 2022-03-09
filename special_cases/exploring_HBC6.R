# This script is meant to investigate what happens with HBC_Chen_5, which has a extremely low sdY, causing adjusted BETA and SE to be extremely big.
library(data.table)

# We applied Undo_scales.R to it, to remove adjusted beta and SE, and to rename original columns accordingly.
ds <- fread("HBC_Chen_32888493_6-hg38.tsv.gz")
n  <- 801305
oneover <- 1/ds$SE^2
nvx <- 2 * n * ds$ALT_FREQ * (1-ds$ALT_FREQ)
m <- lm(nvx ~ oneover - 1)
plot(nvx~oneover - 1)
cf <- coef(m)[['oneover']]

# Sampling for plot
sam10k <- sample(1:length(oneover), 10000, replace=F)
ones <- oneover[sam10k]
nvxs <- nvx[sam10k]
m.s <- lm(nvxs ~ ones - 1)
plot(nvxs~ones-1)
abline(m)

# There's a pesky outlier, what if we got rid of it?
dsnoout <- ds[which(oneover < 1e10), ]
oneover.noout <- 1/dsnoout$SE^2
nvx.noout  <- 2 * n * dsnoout$ALT_FREQ * (1 - dsnoout$ALT_FREQ)
m.noout <- lm(nvx.noout ~ oneover.noout - 1)
plot(nvx.noout ~ oneover.noout - 1)


# Let's compare it to the other meta-analysis (HBC_Chen_..._5)
ds <- fread("../02-Processed/HBC_Chen_32888493_5-hg38.tsv.gz")
ds[,BETA:=orig.BETA][,SE:=orig.SE]

n  <- 801305
oneover <- 1/ds$SE^2
nvx <- 2 * n * ds$ALT_FREQ * (1-ds$ALT_FREQ)
m <- lm(nvx ~ oneover - 1)
cf <- coef(m)[['oneover']]


# Conclusion: HBC_6 has some SNPs with extremely low SE, which in turn makes oneover extremely big, making nvx ~ oneover - 1 slope extremely small. This, in turn makes adjusted BETA and SE very big and the projections get so big that they eclipse other signals.


