###################################
# Computing ALT_FREQ from OConnor #
###################################

# Author: Guillermo Reales
# Date: 2022-03-09

# Background: O'Connor (31189108) comes without allele frequencies, but with genotype counts (n_AA, n_AB, n_BB), so we can compute the frequencies.



library(data.table)
setDTthreads(10)

path <- '~/rds/rds-cew54-basis/02-Processed'
destpath <- '~/rds/rds-cew54-basis/01-Workshop' # We check and copy them back to 02-Processed

fnames  <- dir(path, "OConnor")

for (i in seq_along(fnames)){
print(fnames[i])
dd  <- fread(paste0(path, "/", fnames[i]))
counts <- c("n_AA", "n_AB", "n_BB")
dd[, N:= rowSums(.SD), .SDcols=counts] # Two birds one stone
dd[, ALT_FREQ:= (n_AB/2 + n_BB)/N]     # Allele B is ALT
fwrite(dd, paste0(destpath,"/", fnames[i]), sep = "\t")
}
