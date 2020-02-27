#########################################
## DATASET PROJECTION ONTO THE IMD BASIS
#########################################

# Introduction: This code is meant to project our processed datasets (found at 04-Liftovered/) onto the IMDbasis

##############################################
### LOAD LIBRARIES AND SET REQUIRED FUNCTIONS
##############################################


library(cupcake)
library(data.table)
library(magrittr)

# Change to the right directory

#setwd("~/rds/rds-cew54-wallace-share/Data/Tidy_GWAS-summary/Projecting_dbs/")



#######################
# Load essential files
#######################

# Since files are in hg38 and the basis is in hg19, we use a dirty fix to liftover back to hg19 (easier to do at this stage than in the previous one, due to the huge size of files).
build_dict <- fread("Manifest_build_translator.tsv")
build_dict[, `:=`(pid38=paste(CHR38,BP38,sep=":"), pid=paste(CHR19,BP19,sep=":"))]
build_dict <- build_dict[,c("pid38", "pid")] 

# We load the traits that were used to create the basis first, then we'll append the rest of traits to that table
projected.table <- copy(cupcake::basis.trait.proj)
projected.table[, `:=` (Var.Delta=0, z=0, P=NA)]
setnames(projected.table, c("delta", "trait"), c("Delta", "Trait"))


nfiles <- length(dir("Filtered_datasets/", pattern = ".tsv"))
Trait <- rep(NA,nfiles)
nSNP <- rep(NA, nfiles)
overall_p <- rep(NA, nfiles)
mscomp <- rep(NA,nfiles)


for (i in 1:nfiles){
  
  file <- dir("Filtered_datasets/", pattern = ".tsv")[i]
  trait_label <- strsplit(file, "-", fixed = TRUE)[[1]][1]
  cat("Working on ", file, ".\n", sep = "")
  input_M <- fread(paste("Filtered_datasets/", file, sep = ""))
  input_M <- merge(input_M,build_dict)
  input_M <- na.omit(input_M)

  # A bit of QC
  Trait[i] <- file
  nSNP[i] <- nrow(input_M)
  
  flag <- TRUE
  projected.userdata <- tryCatch(expr = cupcake::project_sparse(beta=input_M$BETA, seb=input_M$SE, pid=input_M$pid)[,trait:=trait_label][],
                        error =function(e) flag <<- FALSE)
  if(!flag){
      message("Projection for the above file had non-zero exit status, please check. Jumping to next file...")
      next 
      }
  projected.userdata[, proj:=NULL]
  setnames(projected.userdata, c("var.proj", "delta", "p", "trait"), c("Var.Delta", "Delta", "P", "Trait")) # Changed var.proj to var.delta
  setcolorder(projected.userdata, c("PC", "Var.Delta", "Delta", "p.overall", "z", "P", "Trait"))
  
  # QC
  overall_p[i] <- sprintf("%1.0e",projected.userdata$p.overall[1])
  minc.idx <- which.min(projected.userdata$P)
  mscomp[i] <- sprintf("%s (%1.0e)",projected.userdata$PC[minc.idx],projected.userdata$P[minc.idx])
  
  

  projected.table <- rbind(projected.table[,.(PC, Delta, Var.Delta, z, P, Trait)],projected.userdata[,.(PC, Delta, Var.Delta, z, P, Trait)])
  
}

QC.table <- data.table(Trait, nSNP, overall_p, mscomp)
QC.table$Trait <- sapply(strsplit(QC.table$Trait, "-"), `[`, 1)

write.table(projected.table, "Projected_table_20200221.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(QC.table, "QC_table_20200221.tsv", sep = "\t", quote = FALSE, row.names = FALSE)



# Tests/debug
#input_M1 <- fread("Filtered_datasets/CD_Liu_26192919_1-ft.tsv")
#input_M1 <- merge(input_M1, build_dict)
#input_M2 <- fread("Filtered_datasets/CD_Liu_26192919_2-ft.tsv")
#input_M2 <- merge(input_M2,build_dict)
#input_M2 <- na.omit(input_M2)
#projected.M1 <- cupcake::project_sparse(beta=input_M1$BETA, seb=input_M1$SE, pid=input_M1$pid)[,trait:="CD1"][]
#projected.M2 <- cupcake::project_sparse(beta=input_M2$BETA, seb=input_M2$SE, pid=input_M2$pid)[,trait:="CD2"][]
#
#input_M <- fread("Filtered_datasets/MDD_Wray_29700475_2-ft.tsv")
#input_M <- merge(input_M, build_dict)
#projected.M <- cupcake::project_sparse(beta=input_M$BETA, seb=input_M$SE, pid=input_M$pid)[,trait:="MDD"][]
#
#
## Project datasets onto the basis
#projected.userdata <- cupcake::project_sparse(beta=input_M$BETA, seb=input_M$SE, pid=input_M$pid)[,trait:="AnkSpond"][]
#projected.userdata[, proj:=NULL] # Removed proj variable
#setnames(projected.userdata, c("var.proj", "delta", "p", "trait"), c("Var.Delta", "Delta", "P", "Trait")) # Changed var.proj to var.delta
#setcolorder(projected.userdata, c("PC", "Var.Delta", "Delta", "p.overall", "z", "P", "Trait"))
#
#basistable.projected <- copy(cupcake::basis.trait.proj)
#setnames(basistable.projected, c("delta", "trait"), c("Delta", "Trait"))
#
#combined.deltaplot.dt <- rbind(basistable.projected[,.(PC,Delta,Var.Delta=0,Trait)],projected.userdata[,.(PC,Delta,Var.Delta,Trait)])
#combined.deltaplot.dt[Var.Delta!=0, ci:=sqrt(Var.Delta) * 1.96]










