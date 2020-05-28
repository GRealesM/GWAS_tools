# This sript is intended to generate Manhattan plots for all files (meta-analyses inputs (2) + meta-analysis)

args <- commandArgs(trailingOnly = TRUE)

library(data.table)
library(ggplot2)
library(magrittr)
library(gridExtra)

trait <- strsplit(args[1], split="_")[[1]][1]
message("Generating plots and table for ", trait)

meta <- fread(args[3])  
input1 <- fread(args[1], select = c("SNPID", "CHR38", "BP38","REF", "ALT", "BETA", "SE", "P"), na.strings="")
input2 <- fread(args[2], select = c("SNPID", "CHR38", "BP38","REF", "ALT", "BETA", "SE", "P"), na.strings="")

# Filter by autosomes
chrs <- as.character(1:22)

meta <- meta[CHR38 %in% chrs,]
input1 <- input1[CHR38 %in% chrs,]
input2 <- input2[CHR38 %in% chrs,]

# Summary table on top 50 SNPs
top20SNPs  <- unique(c(meta[order(P)][1:20,SNPID], input1[order(P)][1:20,SNPID], input2[order(P)][1:20,SNPID])) 
tablemeta  <- meta[SNPID %in% top20SNPs,][,dataset:=args[3]]
tableinput1  <- input1[SNPID %in% top20SNPs,][,dataset:=args[1]]
tableinput2  <- input2[SNPID %in% top20SNPs,][,dataset:=args[2]]
top20table  <- rbind(tablemeta,tableinput1,tableinput2, fill=TRUE)
fwrite <- fwrite(top20table, paste(trait, "_top20table.tsv", sep=""), row.names=FALSE, quote=FALSE, sep = "\t", na = "NA")
message("Table generated!")


manplot <- function(gwas.DT, sig=5e-8, chrs = chrs, title){
	# Plotting so many points is computer-intensive. To save a bit of time and computing resources we'll take all significant SNPs and sample non-significant ones.
	# This way we don't lose much info, while making datasets lighter.
	gwas.sig  <- gwas.DT[P<0.05,]
	gwas.nosig  <- gwas.DT[P>0.05,][sample(.N, .N/5)]
	gwas.DT  <- rbind(gwas.sig,gwas.nosig)
	gwas.DT  <- gwas.DT[order(as.numeric(CHR38),BP38)]
	gwas.DT$CHR38  <- factor(gwas.DT$CHR38, levels=chrs)

	# We need to recodify BP coordenates for chromosomes to appear ordered one after another.
	# BP reference their position in the chromosome. BPcum will inform the "global" coordinates, if BP were counted from chr1 through chr22.
	chr.len  <- gwas.DT[,.(chr.len=max(BP38)),by=CHR38][,tot:=cumsum(as.numeric(chr.len))-chr.len]
	gwas.DT <- merge(gwas.DT,chr.len, by="CHR38")
	gwas.DT[, BPcum:=BP38+tot][, c("chr.len", "tot"):=NULL]
	
	# Some parameters for visualization. 
	# axis.set will place the chr number at the right point in the x axis (middle of the corresponding chromosome)
	# ylim will adjust plot height, leaving some space at the top
	axis.set  <- gwas.DT[,.(center=((min(BPcum)+max(BPcum))/2)), by=CHR38][order(as.numeric(CHR38))]
	ylim <- abs(floor(log10(min(gwas.DT$P)))) + 2

	# Let's prepare the plots
	manhattanplot  <- ggplot(gwas.DT, aes(x = BPcum, y = -log10(P), color = CHR38, size = -log10(P))) +
	  geom_point(alpha = 0.75) +
	  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") +
	  scale_x_continuous(label = axis.set$CHR38, breaks = axis.set$center) +
	  scale_y_continuous(expand = c(0,0), limits = c(0, ylim), breaks=seq(0, ylim, by=2)) +
	  #scale_color_manual(values = rep(c("#276FBF", "#183059"), length(chrs))) +
	  scale_size_continuous(range = c(0.5,3)) +
	  labs(x = NULL,
	       y = "-log10(p)",
	       title = title) +
	  theme_minimal() +
	  theme(legend.position = "none",
		panel.border = element_blank(),
	    panel.grid.major.x = element_blank(),
	    panel.grid.minor.x = element_blank(),
	    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5))
}

metaplot  <- manplot(meta, title=args[3])
inputplot1 <- manplot(input1, title=args[1])
inputplot2  <- manplot(input2, title=args[2])

png(paste("manhattanplot-",trait,".png",sep=""), height=1200, width=800)

grid.arrange(metaplot, inputplot1, inputplot2, nrow=3)

dev.off()
message("Manhattan plots generated!")

