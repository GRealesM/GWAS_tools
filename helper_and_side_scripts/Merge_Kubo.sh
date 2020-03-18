#!/bin/bash
Rscript --vanilla Remerge_Riken.R PTB.auto.rsq07.mac10.txt.gz PTB.chrx.rsq07.mac10.txt.gz TB_Kubo_BBJ_1-raw.tsv;
Rscript --vanilla Remerge_Riken.R NS.auto.rsq07.mac10.txt.gz NS.chrx.rsq07.mac10.txt.gz INS_Kubo_BBJ_1-raw.tsv;
Rscript --vanilla Remerge_Riken.R Pollinosis.auto.rsq07.mac10.txt.gz Pollinosis.chrx.rsq07.mac10.txt.gz ALR_Kubo_BBJ_1-raw.tsv;
Rscript --vanilla Remerge_Riken.R RA.auto.rsq07.mac10.txt.gz RA.chrx.rsq07.mac10.txt.gz RA_Kubo_BBJ_1-raw.tsv;
Rscript --vanilla Remerge_Riken.R T2D.auto.rsq07.mac10.txt.gz T2D.chrx.rsq07.mac10.txt.gz T2D_Kubo_BBJ_1-raw.tsv

