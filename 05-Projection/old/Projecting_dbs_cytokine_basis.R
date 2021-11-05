#########################################
## DATASET PROJECTION ONTO THE BLOOD CELL BASIS
#########################################

# Introduction: This code is meant to project our processed datasets onto the cytokine basis


##############################################
### LOAD LIBRARIES AND SET REQUIRED FUNCTIONS
##############################################

library(data.table)
library(Matrix)
library(magrittr)
load("~/rds/rds-cew54-basis/03-Bases/cytokine_basis/cytokine-basis-sparse-1.0.RData")

# We load the traits that were used to create the basis first, then we'll append the rest of traits to that table
projected.basis <- readRDS("~/rds/rds-cew54-basis/03-Bases/cytokine_basis/cytokine_basis_trait_projection.RDS")
projected.basis[, `:=` (Var.Delta=0, z=0, P=NA)]
setnames(projected.basis, c("delta", "trait"), c("Delta", "Trait"))

files_to_project  <- dir("~/rds/rds-cew54-basis/03-Bases/cytokine_basis/reduced_datasets/", pattern = ".tsv")
nfiles <- length(files_to_project)
Trait <- rep(NA,nfiles)
nSNP <- rep(NA, nfiles)
overall_p <- rep(NA, nfiles)
mscomp <- rep(NA,nfiles)

projected.table <- lapply(files_to_project, function(file){
	message("Projecting ", file)
	trait_label <- strsplit(file, "-", fixed = TRUE)[[1]][1]
	ss.file <- file.path("~/rds/rds-cew54-basis/03-Bases/cytokine_basis/reduced_datasets", file)
	index <- which(files_to_project == file) 
	sm <- fread(ss.file)
	sm <- unique(sm)
	sm <- na.omit(sm)

	# A bit of QC 
	Trait[index] <<- trait_label
	nSNP[index] <<- nrow(sm)

  	projected.userdata <- tryCatch(expr = project_sparse(beta=sm$BETA, seb=sm$SE, pid=sm$pid)[,trait:=trait_label][],
  	                      error =function(e) {
				      message("Projection for the above file had non-zero exit status, please check. Jumping to next file...")
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
projected.full  <- rbind(projected.basis[,.(PC, Delta, Var.Delta, z, P, Trait)],projected.table[,.(PC, Delta, Var.Delta, z, P, Trait)])

QC.table <- data.table(Trait, nSNP, overall_p, mscomp)

date <- format(Sys.time(), format="%Y%m%d")
version  <- 1
projtablename  <- paste("~/rds/rds-cew54-basis/03-Bases/Projections/Projection_cytokine_basis_", date, "-v",version, ".tsv", sep="")
qctablename  <- paste("~/rds/rds-cew54-basis/03-Bases/Projections/QC_cytokine_basis_", date, "-v",version, ".tsv", sep="")

while(file.exists(projtablename)){
  version  <- version + 1
  projtablename  <- paste("~/rds/rds-cew54-basis/03-Bases/Projections/Projection_cytokine_basis_", date, "-v",version, ".tsv", sep="")
  qctablename  <- paste("~/rds/rds-cew54-basis/03-Bases/Projections/QC_cytokine_basis_", date, "-v",version, ".tsv", sep="")
}

write.table(projected.full, projtablename, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(QC.table, qctablename, sep = "\t", quote = FALSE, row.names = FALSE)
