#########################################
## DATASET PROJECTION ONTO THE BLOOD CELL BASIS
#########################################

# Introduction: This code is meant to project our processed datasets onto the Blood cell basis

# We'll use the same reduced datasets as cell basis v2 for testing purposes.


###############################################
######   NEWS                            ######
###############################################

# - On 2021-11-04 I included some lines to identify duplicated pids and remove them, in case they appear

##############################################
### LOAD LIBRARIES AND SET REQUIRED FUNCTIONS
##############################################

library(data.table)
library(Matrix)
library(magrittr)

date <- format(Sys.time(), format="%Y%m%d")
mpath <- "~/rds/rds-cew54-basis/03-Bases/"

load(paste0(mpath, "BCB4/basis_data/BCB4-Astle_1e4_p1e8.RData"))

#Create log file
logname <- paste0(mpath, "BCB4/projections/logs/log_CB4_",date,".txt")
file.create(logname)


files_to_project  <- dir(paste0(mpath,"BCB4/reduced_datasets/"), pattern = ".tsv")
nfiles <- length(files_to_project)
Trait <- rep(NA, nfiles)
nSNP <- rep(NA, nfiles)
overall_p <- rep(NA, nfiles)
mscomp <- rep(NA, nfiles)

projected.table <- lapply(files_to_project, function(file){
	message("Projecting ", file)
	trait_label <- strsplit(file, "-", fixed = TRUE)[[1]][1]
	ss.file <- file.path(paste0(mpath, "BCB4/reduced_datasets"), file)
	index <- which(files_to_project == file) 
	sm <- fread(ss.file) 
	# Some checks
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
	
	ninf <- nrow(sm[BETA == Inf | BETA == -Inf | SE == Inf | SE == -Inf])
    if(ninf > 0) {
		dupmessage= paste0("This file has ", ninf," SNPs with infinite Betas and/or SE. I removed them prior to projection. You might want to check them.")
		message(dupmessage)
		# Write to log
		write(paste0(ss.file,". ",dupmessage), logname, append=TRUE)
		sm <- sm[BETA != Inf & BETA != -Inf & SE != Inf & SE != -Inf] # Remove all infinite values	
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
  	overall_p[index] <<- projected.userdata$p.overall[1] # Note: I updated it to contain the actual number, rather than a reformated one!
  	minc.idx <- which.min(projected.userdata$P)
  	mscomp[index] <<- sprintf("%s (%1.0e)",projected.userdata$PC[minc.idx],projected.userdata$P[minc.idx])
	}
	projected.userdata
}
) 

projected.table <- rbindlist(projected.table, fill = TRUE)
projected.table  <- projected.table[,.(PC, Delta, Var.Delta, z, P, Trait)]

QC.table <- data.table(Trait, nSNP, overall_p, mscomp)

version  <- 1
projtablename  <- paste0("Projection_BCB4_", date, "-v",version, ".tsv")
qctablename  <- paste0("QC_BCB4_", date, "-v",version, ".tsv")

while(projtablename %in% dir(paste0(mpath, "BCB4/projections"))){
  version  <- version + 1
  projtablename  <- paste0("Projection_BCB4_", date, "-v",version, ".tsv")
  qctablename  <- paste0("QC_BCB4_", date, "-v",version, ".tsv")
}

write.table(projected.table, paste0(mpath,"BCB4/projections/", projtablename), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(QC.table, paste0(mpath, "BCB4/projections/", qctablename), sep = "\t", quote = FALSE, row.names = FALSE)

