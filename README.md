# GWAS_tools

In this repository I put all useful scripts that I wrote for processing and projecting GWAS summary statistic datasets onto our multidimensional bases.

Among them, I created a pipeline for initial processing of raw summary statistic datasets. Which changes among versions I log here (from v4 onwards). 


**v4.1** 

* I introduced a 's///' to remove carriage returns from files generated in Windows
* Pipeline now joins liftovered files by keeping only SNPs that succesfully liftovered, as the previous behaviour (keeping all SNPs) was introducing way too many mismatches, including for SNPs that did in fact match. This may result in output files with less lines than input files.

**v4.4**

* Dictionary improved with more terms.
* Pipeline now checks if BETA has more than 50% of SNPs as "NA", and recalculates from OR if it exists and has <50% NA.

**v4.5**

* I included a step in the final step (join, recompress, and save) to check and remove duplicated lines in the file. Because of the join command, some lines were systematically duplicated. This didn't have a huge effect on subsequent steps, but added innecessary lines.

**v4.6**

* SNPID is necessary for liftover step, but it's not absolutely required for it to contain the corresponding rs, so instead of throwing an error and jumping to next file when SNPID column is missing, it creates a new one from CHR and BP columns, with the CHR_BP format. 
* Temporarily included A1 => REF and A2 => ALT to process Hoglund files.
* Added new keys to the dictionary corresponding to Hoglund keys (eg. effB for BETA and se_effB for SE).

**v4.7**

* Fixed a bug that prevented cleaning-up of files when file is in hg38 build.

**v4.8**

* Included some terms in the dictionary to correctly process COVID-19 HGI datasets.

**v4.9**

* Included some terms in the dictionary to correctly process COVID-19 Rivas datasets.
* Removed A1 and A2, included for specific purposes, from the dictionary.
* Now once original build is identified, column names change to reflect it (eg. CHR/BP -> CHR19/BP19), rather than leave them untouched, as it did before.

**v4.10**

* New way to extract file base names, since now some of our datasets include dots in the trait name.

**v5.0**

* New post-liftover merging method, using awk instead of join. This result in seamless merge, keeping all original rows, using NA for missing hg38 coordinates, but without loss of rows due to faulty match, as it happened with the previous version.

