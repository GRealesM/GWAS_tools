# UKBB (Neale) Re-download

**Guillermo Reales
18/07/2022**

We have a new, common manifest to multiple bases, so we'll use the same list of available phenotypes but filtered by the new ~11k SNP manifest (coords_hg19).
Then we'll process them again (ie. liftover, etc.), and include them in our collection.

Note that the files are supposed to be the same, and we're calling this "Neale2" to differentiate the procedure (not the files!) from the first download.
The newly downloaded files will replace the old Neale ones in the collection.

First we'll download the files one by one and filter them by coords_hg19. Note that Neale's coordinate format is a bit strange (header is variant, and format is CHR:BP:REF:ALT) so we adjusted the coordinates accordingly.
Then, we applied the last version of the pipeline (5.3.2_beta) to the files.
This process is included in `slurm_procNeale2`.

Finally, Neale et al. used linear models to generate their GWAS, which is not optimal for binary (ie. case-control) traits (see [here](http://www.nealelab.is/blog/2017/9/11/details-and-considerations-of-the-uk-biobank-gwas)), so we transformed beta and SE from linear to OR scale, using `Fix_scales.R`. This was done by `slurm_fixscaleNeale2`.

