#########################################
## DATASET PROJECTION ONTO THE BLOOD CELL BASIS
#########################################

# Introduction: This code is meant to project our processed datasets onto the Blood cell basis


##############################################
### LOAD LIBRARIES AND SET REQUIRED FUNCTIONS
##############################################

library(data.table)
library(Matrix)
library(magrittr)

date <- format(Sys.time(), format="%Y%m%d")
path = "~/rds/rds-cew54-basis/03-Bases/cell_basis_v2/"
load(paste0(path, "cell-basis-sparse-2.0.RData"))

#Create log file
logname <- paste0(path,"Projections/logs/log_CB2_",date,".txt")
file.create(logname)

# Obsolete
# We load the traits that were used to create the basis first, then we'll append the rest of traits to that table
#projected.basis <- readRDS("../cell_basis_v2/cell_basis_v2_trait_projection.RDS")
#projected.basis[, `:=` (Var.Delta=0, z=0, P=NA)]
#setnames(projected.basis, c("delta", "trait"), c("Delta", "Trait"))

files_to_project  <- dir(paste0(path, "reduced_datasets/"), pattern = ".tsv")
nfiles <- length(files_to_project)
Trait <- rep(NA,nfiles)
nSNP <- rep(NA, nfiles)
overall_p <- rep(NA, nfiles)
mscomp <- rep(NA,nfiles)

projected.table <- lapply(files_to_project, function(file){
	message("Projecting ", file)
	trait_label <- strsplit(file, "-", fixed = TRUE)[[1]][1]
	ss.file <- file.path( paste0(path, "reduced_datasets"), file)
	index <- which(files_to_project == file) 
	sm <- fread(ss.file)
	sm <- unique(sm)
		sm[sm == ""] <- NA # Some missing data might pass as empty string. This will fix that	
	sm <- na.omit(sm, cols = c("pid", "BETA", "SE", "P"))
	dups <- sm$pid[duplicated(sm$pid)]

	if(length(dups) > 0){
		dupmessage= "This file has duplicated pids. I removed them prior to projection. You might want to check it."
		message(dupmessage)
		# Write to log
		write(paste0(ss.file,". ",dupmessage), logname, append=TRUE)
		sm <- sm[!pid %in% dups] # Remove all duplicated instances, to be safe
	}

	# A bit of QC 
	Trait[index] <<- trait_label
	nSNP[index] <<- nrow(sm)

  	projected.userdata <- tryCatch(expr = project_sparse(beta=sm$BETA, seb=sm$SE, pid=sm$pid)[,trait:=trait_label][],
  	                      error =function(e) {
				      failmessage <- "Projection for this file had non-zero exit status, please check. Jumping to next file..."
				      message(failmessage)
				      write(paste0(ss.file,". ",failmessage), logname, append=TRUE)
				      return(NULL)})
  	if(is.null(projected.userdata)) {
		return(NULL)
	} else{
	projected.userdata[, proj:=NULL]
  	setnames(projected.userdata, c("var.proj", "delta", "p", "trait"), c("Var.Delta", "Delta", "P", "Trait")) # Changed var.proj to var.delta
  	setcolorder(projected.userdata, c("PC", "Var.Delta", "Delta", "p.overall", "z", "P", "Trait"))

  	# More QC
  	overall_p[index] <<- projected.userdata$p.overall[1]
  	minc.idx <- which.min(projected.userdata$P)
  	mscomp[index] <<- sprintf("%s (%1.0e)",projected.userdata$PC[minc.idx],projected.userdata$P[minc.idx])
	}
	projected.userdata
}
) 
projected.table[sapply(projected.table, is.null)]  <- NULL 
projected.table <- rbindlist(projected.table)
projected.table  <- projected.table[,.(PC, Delta, Var.Delta, z, P, Trait)]

QC.table <- data.table(Trait, nSNP, overall_p, mscomp)

version  <- 1
projtablename  <- paste("Projection_cell_basis_v2_", date, "-v",version, ".tsv", sep="")
qctablename  <- paste("QC_cell_basis_v2_", date, "-v",version, ".tsv", sep="")

while(projtablename %in% dir(paste0(path, "Projections"))){
  version  <- version + 1
  projtablename  <- paste("Projection_cell_basis_v2_", date, "-v",version, ".tsv", sep="")
  qctablename  <- paste("QC_cell_basis_v2_", date, "-v",version, ".tsv", sep="")
}

write.table(projected.table, paste0(path, "Projections/", projtablename), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(QC.table, paste0(path, "Projections/",qctablename), sep = "\t", quote = FALSE, row.names = FALSE)
