# This script is intended to check variance of Beta in reduced datasets
# Here we'll use metadata to select case-control datasets among reduced datasets for IMD basis,
# then compute N0N1Calc (ie. (N0+N1)/(N0*N1), and plot its log against log(SE) across several SNPs.
# This will allow us to easily identify SNPs and datasets which SE is not proportional to its sample size, and apply further QC to them.
# I selected some SNPs in IMD basis that I know are well known and suited for this (and in different ld.block).

# Guillermo Reales
# 2020-07-30
# Version 1.0

##############################################
### LOAD LIBRARIES AND SET REQUIRED FUNCTIONS
##############################################

library(data.table)
library(googlesheets4)
library(ggplot2)
library(cowplot)
library(ggrepel)

# Load metadata
gs4_deauth()
metadata <- read_sheet("https://docs.google.com/spreadsheets/d/16B4ANehcS4psdAFReTBQUJLYXuf5_RrpmjSPaASa2Nw/edit?usp=sharing", sheet = 1)
metadata <- metadata[,c("File_ID", "Trait_long", "Trait_class", "N0", "N1", "N", "SNP_number", "Chip", "Population", "Collection")]
metadata <- metadata[!is.na(metadata$File_ID),]
metadata$Trait <- sapply(strsplit(metadata$File_ID,"-"),`[`, 1)
metadata <- data.table(metadata)
date <- format(Sys.time(), format="%Y%m%d")

fwrite(metadata,paste0("metadata",date,".tsv"),sep="\t")
#metadata <- fread("metadata20200730.tsv")

reduced.datasets <- sapply(strsplit(dir("~/rds/rds-cew54-basis/03-Bases/IMD_basis/reduced_datasets/"), "-"), `[`,1)

metadata <- metadata[Trait %in% reduced.datasets,]

CC.datasets <- metadata[N1 > 0,]$Trait
CC.metadata <- metadata[Trait %in% CC.datasets, c("Trait","N0","N1")]


reduced_df  <- lapply(CC.datasets, function(x){
	ds  <- fread(paste("~/rds/rds-cew54-basis/03-Bases/IMD_basis/reduced_datasets/", x, "-ft.tsv", sep=""))
	ds$filename  <- x
	ds
}) 
reduced_df  <- rbindlist(reduced_df, fill = TRUE)
reduced_df[,c("alleles","alleles.manifest"):=NULL]
reduced_df <- merge(reduced_df, CC.metadata, by.x="filename", by.y="Trait")

reduced_df$Z <- reduced_df$BETA/reduced_df$SE
reduced_df$N0N1Calc <- (reduced_df$N0 + reduced_df$N1)/(as.numeric(reduced_df$N0) * as.numeric(reduced_df$N1))

SNPs  <- c("rs2476601","rs7130588","rs8067378", "rs11065987", "rs12924729")

# rs2476601 is located in the PTPN22 gene and also known as R620W, C1858T, or 1858C>T, may influence risk for multiple autoimmune diseases, such as Rheumatoid Arthritis, Type-1 diabetes, Autoimmune thyroiditis, and Systemic lupus erythematosus. 
# rs8067378 is associated with childhood asthma
# rs7130588 is a regulatory region variant, associated with asthma and atopic dermatitis.
# rs11065987 is an intergenic variant, associated with platelet count, total cholesterol, and others. [Ensembl](http://www.ensembl.org/Homo_sapiens/Variation/Phenotype?db=core;r=12:111634120-111635120;v=rs11065987;vdb=variation;vf=98111138), [SNPedia](https://www.snpedia.com/index.php/Rs11065987).
# rs12924729** is an intron variant, associated with risk to PBC (G allele). [Ensembl](http://www.ensembl.org/Homo_sapiens/Variation/Phenotype?db=core;r=16:11093426-11094426;v=rs12924729;vdb=variation;vf=24706849).

reduced_df <- reduced_df[SNPID %in% SNPs,]

seplots <- lapply(SNPs, function(snp){
	snp.subset <- reduced_df[SNPID == snp,]
	snp.subset <- na.omit(snp.subset, cols= c("SE","N0N1Calc"))
	snp.subset <- snp.subset[SE != 0,]
	lm.snp.subset <- lm(log(SE) ~ log(N0N1Calc), data = snp.subset)
	cooksd <- cooks.distance(lm.snp.subset)
	sample_size  <- nrow(snp.subset)
	is.outlier <- cooksd > 0.01 # since 4/sample_size was too low 
	
	p  <- ggplot(snp.subset, aes(x = log(N0N1Calc), y = log(SE), colour = BETA, label = filename)) +
	    geom_point()+
	    geom_point(data = snp.subset[is.outlier,], color = "red")+
	    geom_text_repel(size = 2, colour="black", data = snp.subset[is.outlier,], force = 7, box.padding=0.3, seed=4)+
	    labs(title = snp)+
	    theme_classic()
})


zplots <- lapply(SNPs, function(snp){
	snp.subset <- reduced_df[SNPID == snp,]

	z <- ggplot(snp.subset, aes(x = BETA, y = abs(Z), label = filename)) +
	  geom_point()+
	  geom_text_repel(size = 2, colour="black", data = subset(snp.subset,abs(BETA)>0.5*max(abs(BETA))), force = 7, box.padding=0.3, seed=4)+
	  geom_vline(xintercept = 0, linetype = "dashed", colour = "red")+
	  geom_hline(aes(yintercept =5), linetype = "dashed", colour = "red")+
	  labs(title = snp)+
	  theme_classic()
})


plot_grid(plotlist=seplots, ncol=2)
ggsave(paste0("~/rds/rds-cew54-basis/03-Bases/QCplots/SEplots-",date,".png"),height=300,width=200,units="mm")

plot_grid(plotlist=zplots, ncol=2)
ggsave(paste0("~/rds/rds-cew54-basis/03-Bases/QCplots/Zplots-",date,".png"),height=210,width=140,units="mm")


