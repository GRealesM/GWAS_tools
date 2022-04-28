# Comparing Liftovered Berisa structure to updated MacDonald hg38 LD blocks

# Guillermo Reales
# Date: 2022-04-27

# Original Berisa & Pickrell LD blocks (Berisa & Pickrell, 2016) were defined using 1000GP in hg19 build. We tried to simply liftover them and use in basis (and RapidoPGS) work, but this isn't 100% satisfactory. Having LD blocks defined in a panel in hg38 would be ideal, and MacDonald et al., (2022, bioRxiv) did just that using Icelandic individuals. Here we want to compare both distributions of blocks to see how well the old liftovered block distribution approximated the newly computed one.

# Load libraries
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(karyoploteR)

old <- fread("~/rds/rds-cew54-wallace-share/Data/reference/lddetect/old_liftovered/EUR_hg38.bed", header=F)
new <- fread("https://raw.githubusercontent.com/jmacdon/LDblocks_GRCh38/master/data/EUR_LD_blocks.bed")
names(old)  <- names(new)

# Let's see some summaries
summ <-  merge(old[ , .N, by=chr], new[ , .N, by=chr], by="chr", suffixes = c(".old",".new"))
summ[, diff:=N.old-N.new]
# So it seems like the old, liftovered version, defined way more blocks. It's likely that restrcuturing in hg38 has merged some LD blocks together.

# To plot them we need to convert them into GRanges object
og  <- GRanges(seqnames=old$chr, ranges=IRanges(start=old$start, end=old$stop-1), strand="*")
ng  <- GRanges(seqnames=new$chr, ranges=IRanges(start=new$start, end=new$stop-1), strand="*")



pp <- getDefaultPlotParams(plot.type=2)

# Save plot
pdf(file="LiftOver_Berisa_MacDonald_LDblock_comparison.PDF", height=9, width=12)
kpo <- plotKaryotype(genome="hg38", plot.type=2)
kpPlotRegions(kpo, data=og, cex = 0.8, col="deepskyblue", border="grey11")
kpPlotRegions(kpo, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
kpAddLabels(kpo, "Berisa", cex=0.8, col="#888888")
kpAddLabels(kpo, "MacDonald", data.panel=2, cex=0.8, col="#888888")
dev.off()


# Plotting original Berisa (hg19) to see general structure
hg19 <- fread("~/rds/rds-cew54-wallace-share/Data/reference/lddetect/EUR/fourier_ls-all.bed")
hg19 <- GRanges(seqnames=hg19$chr, ranges=IRanges(start=hg19$start, end=hg19$stop-1), strand="*") 
kp19 <- plotKaryotype(genome="hg19")
kpPlotRegions(kp19, data=hg19, cex = 0.8, col="deepskyblue", border="grey11")


# Let's grab a couple of GWAS summary statistic datasets to check if the GWAS peaks fall within the same block defined, which should tell us that the new LD block distribution is ok.
a1 <- fread("~/rds/rds-cew54-basis/02-Processed/AST_Demenais_29273806_1-hg38.tsv.gz")
a1 <- a1[, .(SNPID, CHR38, BP38, P)]
a1[, pid:=paste(CHR38, BP38, sep=":")]
a1 <- a1[!duplicated(pid)]

a1r <- GRanges(seqnames=paste0("chr", a1$CHR38), ranges=IRanges(start=a1$BP38, end=a1$BP38, names=a1$SNPID), strand="*")
a1r$pval  <- a1$P

# This is the full Manhattan plot. We see peaks in multiple chromosomes (eg. 1, 5, 6, 9, and 17). Let's check them one by one.
kpa1 <- plotKaryotype(genome="hg38", plot.type=4)
kpa1 <- kpPlotManhattan(kpa1, data=a1r,  points.col = "brewer.set1")


kpa1 <- plotKaryotype(genome="hg38",	chromosomes=c("chr2","chr5","chr6", "chr9", "chr17"), plot.type=2)
kpa1 <- kpPlotManhattan(kpa1, data=a1r,  points.col = "brewer.set1")
kpPlotRegions(kpa1, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
kpAddLabels(kpa1, "MacDonald", data.panel=2, cex=0.8, col="#888888")

# Let's look closer at some other chromosomes
kpa1 <- plotKaryotype(genome="hg38",	chromosomes=c("chr8","chr10","chr11", "chr15", "chr16"), plot.type=2)
kpa1 <- kpPlotManhattan(kpa1, data=a1r,  points.col = "brewer.set1")
kpPlotRegions(kpa1, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
kpAddLabels(kpa1, "MacDonald", data.panel=2, cex=0.8, col="#888888")


# There's a peak in Chr6 outside of HLA that seems to be between two blocks, let's find out what is it
a1[CHR38 == 6 & P < 5e-8][order(BP38)]
# It looks like the ~6:90000000 region.
kpa1f1 <- plotKaryotype(genome="hg38",	zoom="chr6:80e6-91e6", plot.type=2)
kpa1f1 <- kpPlotManhattan(kpa1f1, data=a1r,  points.col = "brewer.set1")
kpAddBaseNumbers(kpa1f1, add.units = TRUE, tick.dist=1000000, cex=1)
kpPlotRegions(kpa1f1, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
kpAddLabels(kpa1f1, "MacDonald", data.panel=2, cex=0.8, col="#888888")

# There's a peak in Chr10 that seems to be between two blocks, let's find out what is it
a1[CHR38 == 10 & P < 5e-8][order(BP38)]
# It looks like the ~10:9000000 region.
kpa1f2 <- plotKaryotype(genome="hg38",	zoom="chr10:8.8e6-9.2e6", plot.type=2)
kpa1f2 <- kpPlotManhattan(kpa1f2, data=a1r,  points.col = "brewer.set1")
kpAddBaseNumbers(kpa1f2, add.units = TRUE, tick.dist=1000000, cex=1)
kpPlotRegions(kpa1f2, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
kpAddLabels(kpa1f2, "MacDonald", data.panel=2, cex=0.8, col="#888888")

# There's a peak in Chr5 that seems to be between two blocks, let's find out what is it
a1[CHR38 == 5 & P < 5e-8][order(BP38)]
# It looks like the ~5:110M-114M region.
kpa1f3 <- plotKaryotype(genome="hg38",	zoom="chr5:110e6-138e6", plot.type=2)
kpa1f3 <- kpPlotManhattan(kpa1f3, data=a1r,  points.col = "brewer.set1")
kpAddBaseNumbers(kpa1f3, add.units = TRUE, cex=1)
kpPlotRegions(kpa1f3, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
kpAddLabels(kpa1f3, "MacDonald", data.panel=2, cex=0.8, col="#888888")

kpa1f3 <- plotKaryotype(genome="hg38",	zoom="chr5:130e6-135e6", plot.type=2)
kpa1f3 <- kpPlotManhattan(kpa1f3, data=a1r,  points.col = "brewer.set1")
kpAddBaseNumbers(kpa1f3, add.units = TRUE, cex=1)
#kpPlotRegions(kpa1f3, data=og, cex = 0.8, col="deepskyblue", data.panel=2, border="grey11")
kpPlotRegions(kpa1f3, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
#kpAddLabels(kpa1f3, "Berisa", cex=0.8, data.panel=2, col="#888888")
kpAddLabels(kpa1f3, "MacDonald", data.panel=2, cex=0.8, col="#888888")


# Let's try a different dataset
ra1 <- fread("~/rds/rds-cew54-basis/02-Processed/RA_Okada_24390342_1-hg38.tsv.gz")
ra1 <- ra1[, .(SNPID, CHR38, BP38, P)]
ra1[, pid:=paste(CHR38, BP38, sep=":")]
ra1 <- ra1[!duplicated(pid)]
ra1 <- na.omit(ra1)

ra1r <- GRanges(seqnames=paste0("chr", ra1$CHR38), ranges=IRanges(start=ra1$BP38, end=ra1$BP38, names=ra1$SNPID), strand="*")
ra1r$pval  <- ra1$P

# This is the full Manhattan plot. We see peaks in multiple chromosomes (eg. 1, 5, 6, 9, and 17). Let's check them one by one.
kpra1 <- plotKaryotype(genome="hg38", plot.type=4)
kpra1 <- kpPlotManhattan(kpra1, data=ra1r,  points.col = "brewer.set3", ymax=28)
kpAxis(kpra1, ymin = 0, ymax=28)

# Change chr to see different regions
kpra1 <- plotKaryotype(genome="hg38",	chromosomes=c("chr22"), plot.type=2)
kpra1 <- kpPlotManhattan(kpra1, data=ra1r,  points.col = "brewer.set1", ymax=20)
kpAxis(kpra1, ymin = 0, ymax=20)
kpPlotRegions(kpra1, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
kpAddBaseNumbers(kpra1, add.units = TRUE, cex=1)
kpAddLabels(kpra1, "MacDonald", data.panel=2, cex=0.8, col="#888888")

# Let's zoom 
# Suspicious region in chr1
kpra1f1 <- plotKaryotype(genome="hg38",	zoom="chr1:110e6-117e6", plot.type=2)
kpra1f1 <- kpPlotManhattan(kpra1f1, data=ra1r,  points.col = "brewer.set1")
kpAddBaseNumbers(kpra1f1, add.units = TRUE, cex=1)
kpPlotRegions(kpra1f1, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
kpAddLabels(kpra1f1, "MacDonald", data.panel=2, cex=0.8, col="#888888")

# Another suspicious region in chr4
kpra1f1 <- plotKaryotype(genome="hg38",	zoom="chr4:22e6-30e6", plot.type=2)
kpra1f1 <- kpPlotManhattan(kpra1f1, data=ra1r,  points.col = "brewer.set1", ymax=20)
kpAxis(kpra1, ymin = 0, ymax=20)
kpAddBaseNumbers(kpra1f1, add.units = TRUE, cex=1)
kpPlotRegions(kpra1f1, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
kpAddLabels(kpra1f1, "MacDonald", data.panel=2, cex=0.8, col="#888888")

# Another suspicious region in chr5
kpra1f1 <- plotKaryotype(genome="hg38",	zoom="chr5:55e6-60e6", plot.type=2)
kpra1f1 <- kpPlotManhattan(kpra1f1, data=ra1r,  points.col = "brewer.set1", ymax=20)
kpAxis(kpra1, ymin = 0, ymax=20)
kpAddBaseNumbers(kpra1f1, add.units = TRUE, cex=1)
kpPlotRegions(kpra1f1, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
kpAddLabels(kpra1f1, "MacDonald", data.panel=2, cex=0.8, col="#888888")


# Another suspicious region in chr10
kpra1f1 <- plotKaryotype(genome="hg38",	zoom="chr10:5e6-10e6", plot.type=2)
kpra1f1 <- kpPlotManhattan(kpra1f1, data=ra1r,  points.col = "brewer.set3", ymax=20)
kpAxis(kpra1, ymin = 0, ymax=20)
kpAddBaseNumbers(kpra1f1, add.units = TRUE, cex=1)
kpPlotRegions(kpra1f1, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
kpAddLabels(kpra1f1, "MacDonald", data.panel=2, cex=0.8, col="#888888")

# Another suspicious region in chr19
kpra1f1 <- plotKaryotype(genome="hg38",	zoom="chr19:18e6-21e6", plot.type=2)
kpra1f1 <- kpPlotManhattan(kpra1f1, data=ra1r,  points.col = "brewer.set1", ymax=20)
kpAxis(kpra1, ymin = 0, ymax=20)
kpAddBaseNumbers(kpra1f1, add.units = TRUE, cex=1)
kpPlotRegions(kpra1f1, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
kpAddLabels(kpra1f1, "MacDonald", data.panel=2, cex=0.8, col="#888888")

# Another suspicious region in chr20
kpra1f1 <- plotKaryotype(genome="hg38",	zoom="chr20:45e6-48e6", plot.type=2)
kpra1f1 <- kpPlotManhattan(kpra1f1, data=ra1r,  points.col = "brewer.set1", ymax=20)
kpAxis(kpra1, ymin = 0, ymax=20)
kpAddBaseNumbers(kpra1f1, add.units = TRUE, cex=1)
kpPlotRegions(kpra1f1, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
kpAddLabels(kpra1f1, "MacDonald", data.panel=2, cex=0.8, col="#888888")


# Another suspicious region in chr22
kpra1f1 <- plotKaryotype(genome="hg38",	zoom="chr22:35e6-40e6", plot.type=2)
kpra1f1 <- kpPlotManhattan(kpra1f1, data=ra1r,  points.col = "brewer.set1", ymax=20)
kpAxis(kpra1, ymin = 0, ymax=20)
kpAddBaseNumbers(kpra1f1, add.units = TRUE, cex=1)
kpPlotRegions(kpra1f1, data=ng, cex = 0.8, col="gold", data.panel=2, border="grey11")
kpAddLabels(kpra1f1, "MacDonald", data.panel=2, cex=0.8, col="#888888")




