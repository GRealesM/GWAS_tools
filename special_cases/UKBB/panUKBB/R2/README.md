# PanUKBB (Release 2) proceedings

**Guillermo Reales
18/07/2022**

PanUKBB, a project that performed GWASs over 7,000+ phenotypes using UKBB data, was released on 16/06/2020.
On April 2022 a new set of files (which we'll call R2 from now on) were released. Among other things, the authors re-computed the meta-analyses, by applying 'stringent QC' procedures (see [here](https://pan.ukbb.broadinstitute.org/docs/qc)). This resulted in new 'hq' columns, resulting from running the meta-analyses on HQ individuals and variants. 
Like we did with Neale files and the previous PanUKBB release (R1) we intend to download some files of interest to us and filter them by our SNP manifest. This time, we'll use the common [11k SNP manifest](https://github.com/GRealesM/Bases/blob/master/Manifests/Consensus_interbasis_manifest.R).
For consistency, I selected the same phenotypes as those in PanUKBBR1 (1,553 endpoints),  from the [PanUKBB manifest](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=268241601), which resulted 1,552 phenotypes (ph613.9, Breast disorder NOS was unavailable in R2), including 20002 traits, COVID traits, biomarkers, some continuous traits and phecodes.

Since our interest is still to extract as much power as we can, when a meta-analysis is available, we will consider meta-analysis values, rather than individual population values.

The step prior to pipelining is to rehead each file appropriately. For that I extracted the headers for each file header (`extract-headers.sh`), and classified them according to the type of dataset (ie. meta_hq, meta_nohq, EUR, AFR) in `files_by_header`.

Then using this classification, I replaced the right beta, se, and p-value by BETA, SE, P. After choosing and renaming the right estimates for each dataset, I pipelined them (5.3.2_beta). This was done in `slurm_proc_PanUKBBR2`.

