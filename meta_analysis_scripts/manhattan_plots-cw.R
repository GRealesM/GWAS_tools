# This sript is intended to generate Manhattan plots for all files (meta-analyses inputs (2) + meta-analysis)

args <- commandArgs(trailingOnly = TRUE)

library(data.table)
library(ggplot2)
library(magrittr)
library(gridExtra)

#trait <- strsplit(args[1], split="_")[[1]][1]
trait="MIG" # for testing
message("Generating plots and table for ", trait)
args <- c("bNGF_AholaOlli_27989323_1-hg38.tsv.gz",
          "bNGF_Hoglund_31727947_1-hg38.tsv.gz",
          "bNGF_Reales_up_1-hg38.tsv.gz")
setwd("~/share2/04-Cytokine_metaanalysis")
meta <- fread(args[3])  
input1 <- fread(args[1], select = c("SNPID", "CHR38", "BP38","REF", "ALT", "BETA", "SE", "P"), na.strings="")
input2 <- fread(args[2], select = c("SNPID", "CHR38", "BP38","REF", "ALT", "BETA", "SE", "P"), na.strings="")

# Filter by autosomes
chrs <- as.character(1:22)

meta <- meta[CHR38 %in% chrs,]
input1 <- input1[CHR38 %in% chrs,]
input2 <- input2[CHR38 %in% chrs,]

## check on numbers
nrow(input1)
nrow(input2)
nrow(meta)  

## keep snps in all datasets - otherwise no worthwhile comparison
## input1[,PID:=paste(CHR38,BP38,sep=":")]
## input2[,PID:=paste(CHR38,BP38,sep=":")]
## meta[,PID:=paste(CHR38,BP38,sep=":")]
dups <- c(meta$SNPID[ duplicated(meta$SNPID) ],
          input1$SNPID[ duplicated(input1$SNPID) ],
          input2$SNPID[ duplicated(input2$SNPID) ])
snps.keep <- setdiff(meta$SNPID,dups)

meta <- meta[ (SNPID %in% snps.keep) ]
input1 <- input1[SNPID %in% snps.keep]
input2 <- input2[SNPID %in% snps.keep]
## input1 <- input1[SNPID %in% meta$SNPID]
## input2 <- input2[SNPID %in% meta$SNPID]

## ## check again on numbers
## nrow(input1) # 8438997
## nrow(input2) # 8432585
## nrow(meta) # 8438997
## ## Still appears to be some difference

## ## Question: why are these not equal?
## any(duplicated(input1$SNPID)) # TRUE
## any(duplicated(meta$SNPID)) # TRUE
## dups <- c(meta$SNPID[ duplicated(meta$SNPID) ],
##           input1$SNPID[ duplicated(input1$SNPID) ],
##           input2$SNPID[ duplicated(input2$SNPID) ])
## meta <- meta[ !(SNPID %in% dups) ]
## input1 <- input1[SNPID %in% meta$SNPID]
## input2 <- input2[SNPID %in% meta$SNPID]

## check again on numbers
nrow(input1) #8428591
nrow(input2) #8428591
nrow(meta) #8428591

## check effects aligned
check <- input1[P<1e-6]$SNPID
input1[SNPID %in% check]
input2[SNPID %in% check]
meta[SNPID %in% check]

## # Summary table on top 50 SNPs
## top20SNPs  <- unique(c(meta[order(P)][1:20,SNPID], input1[order(P)][1:20,SNPID], input2[order(P)][1:20,SNPID])) 
## tablemeta  <- meta[SNPID %in% top20SNPs,][,dataset:=args[3]]
## tableinput1  <- input1[SNPID %in% top20SNPs,][,dataset:=args[1]]
## tableinput2  <- input2[SNPID %in% top20SNPs,][,dataset:=args[2]]
## top20table  <- rbind(tablemeta,tableinput1,tableinput2, fill=TRUE)
## fwrite <- fwrite(top20table, paste(trait, "_top20table.tsv", sep=""), row.names=FALSE, quote=FALSE, sep = "\t", na = "NA")
## message("Table generated!")

sigsnps <- c(meta[P<1e-2]$SNP, input1[P<1e-2]$SNP, input2[P<1e-2]$SNP)  %>% unique()
nsigsnps <- setdiff(meta$SNP,sigsnps)
nsigsnps %<>% sample(., length(nsigsnps)/100)
usesnps <- c(sigsnps,nsigsnps)

library(cowplot)
manplot <- function(gwas.DT, sig=5e-8, chrs = chrs, title){
	# Plotting so many points is computer-intensive. To save a bit of time and computing resources we'll take all significant SNPIDs and sample non-significant ones.
	# This way we don't lose much info, while making datasets lighter.
	## gwas.sig  <- gwas.DT[P<0.05,]
	## gwas.nosig  <- gwas.DT[P>0.05,][sample(.N, .N/5)]
	gwas.DT  <- gwas.DT[SNPID %in% usesnps]
	gwas.DT  <- gwas.DT[order(as.numeric(CHR38),BP38)]
	## gwas.DT$CHR38  <- factor(gwas.DT$CHR38, levels=chrs)

	# We need to recodify BP coordenates for chromosomes to appear ordered one after another.
	# BP reference their position in the chromosome. BPcum will inform the "global" coordinates, if BP were counted from chr1 through chr22.
	## chr.len  <- gwas.DT[,.(chr.len=max(BP38)),by=CHR38][,tot:=cumsum(as.numeric(chr.len))-chr.len]
	## gwas.DT <- merge(gwas.DT,chr.len, by="CHR38")
	## gwas.DT[, BPcum:=BP38+tot][, c("chr.len", "tot"):=NULL]

  # I think this is simpler than a merge
  gwas.DT[,dpos:=c(0,diff(BP38)),by="CHR38"]
  gwas.DT[dpos<0,dpos:=0]
  gwas.DT[,BPcum:=cumsum(dpos)]
  
	# Some parameters for visualization. 
	# axis.set will place the chr number at the right point in the x axis (middle of the corresponding chromosome)
	# ylim will adjust plot height, leaving some space at the top
	axis.set  <- gwas.DT[,.(center=((min(BPcum)+max(BPcum))/2)), by=CHR38][order(as.numeric(CHR38))]
	ylim <- abs(floor(log10(min(gwas.DT$P)))) + 2

	# Let's prepare the plots
  manhattanplot  <- ggplot(gwas.DT,
                           aes(x = BPcum,
                               y = pmin(-log10(P),20),
                               color = as.numeric(CHR38) %% 2 ==0#, size = -log10(P)
                               )) +
	  geom_point(alpha = 0.75,size=1) +
	  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") +
	  scale_x_continuous(label = axis.set$CHR38, breaks = axis.set$center) +
	  scale_y_continuous(expand = c(0,0), limits = c(0, ylim), breaks=seq(0, ylim, by=2)) +
          scale_colour_grey() +
	  #scale_color_manual(values = rep(c("#276FBF", "#183059"), length(chrs))) +
	  ## scale_size_continuous(range = c(0.5,3)) +
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

plot_grid(metaplot, inputplot1, inputplot2, nrow=3)

dev.off()
message("Manhattan plots generated!")

