# GWAS_Sumstats_Pipeline
Pipeline for GWAS summary statistic dataset processing 

From v4 on I started logging changes across versions, just in case.

**v4.1** 

* I introduced a 's///' to remove carriage returns from files generated in Windows
* Pipeline now joins liftovered files by keeping only SNPs that succesfully liftovered, as the previous behaviour (keeping all SNPs) was introducing way too many mismatches, including for SNPs that did in fact match. This may result in output files with less lines than input files.

**v4.4**

* Dictionary improved with more terms.
* Pipeline now checks if BETA has more than 50% of SNPs as "NA", and recalculates from OR if it exists and has <50% NA.

**v4.5**

* I included a step in the final step (join, recompress, and save) to check and remove duplicated lines in the file. Because of the join command, some lines were systematically duplicated. This didn't have a huge effect on subsequent steps, but added innecessary lines.
