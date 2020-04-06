#!/usr/bin/env Rscript

setwd("~/Projects/ext_basis")
library(data.table)
library(magrittr)
library(parallel)

## sets various file locations
ROOT.DIR <- '~/share/as_basis/basis-creation'
REF_GT_DIR <- file.path(ROOT.DIR,'../GWAS/ctrl_gt/by.chr')

## handy functions
source("cw-utils.R")

## get pids
pids <- scan("~/share/Data/Tidy_GWAS-summary/sparse.cellbasis.SNP.hg19.txt",what="")
chr <- sub(":.*","",pids)
pids.bychr <- split(pids,chr)

## LD
LD <- lapply(names(pids.bychr), function(chr) {
        ss.file <- file.path(REF_GT_DIR, sprintf("%s.RDS", chr))
        sm <- readRDS(ss.file)
        pids <- colnames(sm)
        dup.idx <- which(duplicated(pids))
        if (length(dup.idx) > 0) {
            sm <- sm[, -dup.idx]
            pids <- pids[-dup.idx]
        }
        sm.map <- match(pids.bychr[[chr]], pids)
        if (any(is.na(sm.map))) {
            message("Pids in manifest that don't have genotypes")
        }
        r <- ld(sm[, sm.map], sm[, sm.map], stats = "R")
        r[is.na(r)] <- 0
        r
})  %>% bdiag_with_dimnames(.)

LD <- LD[pids,pids]

## todo: reorder to match rot.pca

## todo: save
