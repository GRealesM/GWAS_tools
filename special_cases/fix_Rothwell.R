## Fix Rothwell datasets 


# Rothwell datasets (Myositis and subtypes) have their chromosomes 1-9 named as 01-09. 
# This produces an error when liftovering, so we need to make them be 1-9.
# The easiest way to do this is to let R read the file and write it as the gods intended.

library(data.table)
setDTthreads(10) # Look at me watchdog!

files = dir(pattern = "*Rothwell_up_1-raw*") # Assuming the datasets are in the same directory

for(i in files){ 
	message(i)
	x = fread(i)
	fwrite(x=x, file=i, sep="\t", na=NA, quote=FALSE) # Keep NA as NA
}

# And that should be it
