# GWAS_Sumstats_Pipeline
Pipeline for GWAS summary statistic dataset processing 

From v4 on I started logging changes across versions, just in case.

**v4.1** 

* I introduced a 's///' to remove carriage returns from files generated in Windows
* Pipeline now joins liftovered files by keeping only SNPs that succesfully liftovered, as the previous behaviour (keeping all SNPs) was introducing way too many mismatches, including for SNPs that did in fact match. This may result in output files with less lines than input files.
