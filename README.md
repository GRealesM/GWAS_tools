# GWAS_tools

This repository is intended to harbour all GWAS summary statistics processing-related code.
GWAS summary statistic datasets can be a bit of a headache to process because they come in many different formats, and builds.
Our goal is to have a set of GWAS summary statistics datasets with a minimum set of columns, including genomic coordinates in hg38,
reference and effect allele, Beta, SE of Beta, and P-values, as well as other informative columns that the file may provide.

We decided to write the following processing and QC protocol, from pre-download to projected dataset onto the bases, for our inner control,
and to allow for reproducibility.

Bear in mind that 
1. This is not an automatic pipeline, so files will need to be carefully checked by the researcher for possible errors and needed fixes.
2. Some of the required files to perform the fixes (eg. reference panels) are too big for Github, so we're not bundling them here.

We structured the code in different steps, and we'll include the code for each section within the correspondingly numbered directories.

**Step 0 - Prior to processing**

Before downloading:

- Check if file is public and downloadable. If not, consider writing to the authors.
- Gather minimum required information about the dataset (Author, PubMedID, Year, Trait, N0, N1, build, array, download link), and add it to table.

After download:

- Check if datasets for the same trait come in multiple file. If so, merge them.
- Check if essential column names are identifiable. This means that all essential columns (see below) must be identifiable by the dictionary in the pipeline, so it's good practice to take a look at the file and see if it has them. If not, pipeline will complain. If some column names are not clear, check original paper or Readme file, if you're lucky enough to have one with your file.
Special note on alleles: REF and ALT alleles are named different ways (eg. A1/A2, A2/A1, a0/a1, AlleleA/AlleleB) which can lead to ambiguitity, so it's a good practice to manually check and replace by REF/ALT if necessary.
- Add new column names to the dictionary in the pipeline (see `01-Main_pipeline/Full_pipeline_vX.X.sh`). Upgrade Pipeline version if so.
- If not a .gz already, compress using `gzip`.

- **0a Missing coordinates** If genomic coordinates are missing, pipeline will fail. If we have SNPIDs, we can use (...) to extract hg19 coordinates from 1000 Genomes project phase III. Ensure your SNP id column is recoded as SNPID. Output will have CHR/BP columns, rather than CHR19/BP19, as this will be done at processing step.
- **0b - Missing one allele** Some authors report only one allele (usually the risk/effect allele), but for many analysis we need both alleles. Use `fix-one-allele-problem.sh` to use 1000 Genomes phase III reference panel (not supplied) to guess the other allele. If you have multiple files with a common list of alleles, you can use `fix-one-allele-problem_commonSNPs.sh`.

- **0c - Misplaced coords** Some (especially old) files report BP as zero-based, so BP is always one base below actual coordinates. This can be fixed using `Update_positions.R`.

- If the file is missing effect size (encoded as OR, BETA, or Z), the file must be discarded, as it doesn't contain enough information for us to proceed.


**Step 1 - Processing file using Pipeline (v5.2)**

`Full_pipeline_vX.X.sh` will perform the following steps automatically, taking a compressed (.tsv.gz) as an input
- Check column separator, in the dataset, and change it to tabs.
- Use dictionary to replace column names in many flavours to the following:
    - Chromosome (CHR)
    - Position (BP)
    - SNP id (SNPID)
    - Reference allele (REF)
    - Alternative, or effect allele (ALT)
    - P-value for association (P)
    - Odds-ratio (OR) or BETA (BETA)
    - Standard Error of the log OR (SE)
    - Effect allele frequency (ALT_FREQ)
- Check if all minimum columns are present. CHR, BP, and P are absolutely essential. OR or BETA must be present, too.
- Check alleles are expressed in upper case letters, and coerce them to upper case if not.
- If columns are missing, try to compute them:
    - Compute BETA from log(OR), if OR is present.
    - Compute BETA and SE from Z, MAF, and N0/N1, if those are present.*
    - Compute SE if missing, if BETA and Z (or P) are present.
    - Create SNPID from CHR_BP, if missing.
    - If only REF or ALT is available, use 1000 Genomes Phase III to compute the other allele. Step not available in pipeline. Use `fix-one-allele-problem.sh`.*
    - For certain steps, we may need MAF. If so, and if MAF is not available in the file, compute them using 1000 Genome Phase III reference and the appropriate proxy population. Use `Compute_freqs.sh`*
- Check missing data for P, BETA or OR column is <50%.
- Replace "23" by "X" at CHR column.
- Use SNP manifest in cupcake basis to identify dataset genome build, by majority rule (hg18, hg19, or hg38).
- If not in hg38, use liftOver executable and appropriate chain files to create two additional columns with hg38 coordinates. Also, rename CHR and BP to CHRxx and BPxx, depending on the build of the original data.

*Steps not implemented in the pipeline yet. See fixes below

**Step 2 - Additional fixes**

Sometimes the files don't have what we need, so we've included some fixes for datasets.
- **2a - Missing allele frequencies** If MAF are required, but not included in the file, `Compute_freqs.sh` can be used to extract allele frequencies from the 1000 Genomes phase III reference panel (not supplied!). Current usage is `Compute_freqs.sh -f file.tsv.gz -p POPULATION`, where POPULATION is the population code for one of the 1000 Genomes Phase III populations. See [here](https://www.internationalgenome.org/faq/which-populations-are-part-your-study/). Note that 1000G reference panel is in hg19 build, so file should include CHR19 and BP19 columns. Future implementations will include taking multiple files in one go.
- **2b - Compute BETA and SE from Z** When BETA and SE is missing in case-control datasets, we can compute them from Z scores by using ALT frequencies and proportion of cases/controls, if available, using (...).
- **2c - Backliftover to hg19** Although we try to use hg38 as gold standard, for many purposes hg19 are necessary. Pipeline can identify build (hg18, hg19, and hg38, currently) and rename CHR/BP columns to their corresponding build. However, if CHR19/BP19 are needed but not, we can *backliftover* to hg19. Use (...).
 

**Step 3 - Further QC**

- **3a Plot SE vs. log(N0+N1/(N0*N1)) for Case-control datasets** Case-control datasets should be in log(OR) scale, and SE should be inversely proportional to sample size. By plotting log(SE) vs. log(N0+N1/(N0*N1)) on several, pre-defined SNPs, they should show as a straight line. Deviations from this line might mean problems with the dataset, like CC being in linear scale (See 3b for solution). In practice, we'll use Reduced files to generate these plots, since those are lighter.

- **3b - Linear to log(OR) scaling in Case-control datasets** CC datasets should have their BETAs and SE in the log(OR) scale, rather than linear scale. Datasets generated using linear mixed models (eg. BOLT-LMM) should be transformed to log(OR) scale prior to procesing. For scale conversion, use `IMDtools::linORscale`.

-**3c Adjust sdY for quantitative datasets** BETA and SE can measure different things across quantitative trait dataset. To make them comparable, they must have variance ≃ 1. Estimate sdY, and divide BETA and SE by it to ensure variance = 1. This can be done using `IMDtools::sdy.correction`. Note that MAF is required for this step.

- A script for performing steps 03b and 03c, depending on the type of dataset is available at `03bc-Convert_scales_adjust_sdY/Fix_scales_adjust_sdY.R`.

- Files processed using the previous 3 steps should be placed in 02-Processed.

**Step 4 - Reduction**

- Prior to projection, we will "reduce" files. This consists in using a basis manifest, filter processed sumstat dataset by the few SNPs in the corresponding basis, align alleles to manifest, and store it in `reduced_datasets/` (within each basis directory in `03-Bases`). Use `Reducing_for_*.R` for that. In practice, we can use those reduced datasets to generate SEvsN0N1Calc (step 3a) plots, as they're way lighter than original files.

**Step 5 - Projection**

- Project files onto the bases. This will use `project_sparse` function within each basis `RData` file (or `cupcake`, if IMD basis), to project each reduced dataset onto the basis, outputting two data.frames: (1) **Projections**, containing Delta, Var.Delta, Z, and P for each PC in the basis, and (2) **QC table**, containing nSNPs matched between the dataset and the basis, overall P-value from the projection, and the lowest-P PC. The outputs will be saved to `03-Bases/Projections` directory.



