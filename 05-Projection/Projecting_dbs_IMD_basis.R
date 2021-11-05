#########################################
## DATASET PROJECTION ONTO THE IMD BASIS
#########################################

# Introduction: This code is meant to project our processed datasets (found at 04-Liftovered/) onto the IMDbasis

###############################################
######   NEWS                            ######
###############################################

# - On 2021-11-04 I included some lines to identify duplicated pids and remove them, in case they appear. I also replace empty strings by NA.

##############################################
### LOAD LIBRARIES AND SET REQUIRED FUNCTIONS
##############################################


library(cupcake)
library(data.table)
library(magrittr)

date <- format(Sys.time(), format="%Y%m%d")
mpath <- "~/rds/rds-cew54-basis/03-Bases/"

#Create log file
logname <- paste0("logs/log_IMD_",date,".txt")
file.create(logname)

#######################
# Load essential files
#######################

# Since files are in hg38 and the basis is in hg19, we use a dirty fix to liftover back to hg19 (easier to do at this stage than in the previous one, due to the huge size of files).
build_dict <- fread(paste0(mpath, "IMD_basis/Manifest_build_translator.tsv"))
build_dict[, `:=`(pid38=paste(CHR38,BP38,sep=":"), pid=paste(CHR19,BP19,sep=":"))]
build_dict <- build_dict[,c("pid38", "pid")] 

# We load the traits that were used to create the basis first, then we'll append the rest of traits to that table
projected.basis <- copy(cupcake::basis.trait.proj)
projected.basis[, `:=` (Var.Delta=0, z=0, P=NA)]
setnames(projected.basis, c("delta", "trait"), c("Delta", "Trait"))

files_to_project  <- dir(paste0(mpath,"IMD_basis/reduced_datasets/"), pattern = ".tsv")
nfiles <- length(files_to_project)
Trait <- rep(NA,nfiles)
nSNP <- rep(NA, nfiles)
overall_p <- rep(NA, nfiles)
mscomp <- rep(NA,nfiles)

projected.table <- lapply(files_to_project, function(file){
	message("Projecting ", file)
	trait_label <- strsplit(file, "-", fixed = TRUE)[[1]][1]
	ss.file <- file.path(paste0(mpath,"IMD_basis/reduced_datasets"), file)
	index <- which(files_to_project == file) 
	sm <- fread(ss.file)
	# Some checks
	sm <- unique(sm)
	sm[sm == ""] <- NA # Some missing data might pass as empty string. This will fix that	
	sm <- na.omit(sm, cols = c("pid38", "BETA", "SE", "P"))
	dups <- sm$pid38[duplicated(sm$pid38)]

	if(length(dups) > 0){
		dupmessage= "This file has duplicated pids. I removed them prior to projection. You might want to check it."
		message(dupmessage)
		# Write to log
		write(paste0(ss.file,". ",dupmessage), logname, append=TRUE)
		sm <- sm[!pid38 %in% dups] # Remove all duplicated instances, to be safe
	}
	
	sm <- merge(sm, build_dict)

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
projected.full  <- rbind(projected.basis[,.(PC, Delta, Var.Delta, z, P, Trait)],projected.table[,.(PC, Delta, Var.Delta, z, P, Trait)])

QC.table <- data.table(Trait, nSNP, overall_p, mscomp)

version  <- 1
projtablename  <- paste0(mpath,"IMD_basis/Projections/Projection_IMD_basis_", date, "-v",version, ".tsv")
qctablename  <- paste0(mpath, "IMD_basis/Projections/QC_IMD_basis_", date, "-v",version, ".tsv")

while(file.exists(projtablename)){
  version  <- version + 1
  projtablename  <- paste0(mpath, "IMD_basis/Projections/Projection_IMD_basis_", date, "-v",version, ".tsv")
  qctablename  <- paste0(mpath, "IMD_basis/Projections/QC_IMD_basis_", date, "-v",version, ".tsv")
}

write.table(projected.full, projtablename, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(QC.table, qctablename, sep = "\t", quote = FALSE, row.names = FALSE)

