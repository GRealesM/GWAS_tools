"%notin%" <- function(x, table) match(x, table, nomatch = 0) == 0

library(data.table)
dt2mat <- function(dt,...) {
    tmp <- dcast(dt,...)
    rn <- tmp[[1]]
    m <- as.matrix(tmp[,-1])
    rownames(m) <- rn
    m
}

indep_snps <- function(dt,LD,thin=TRUE,rthr=0.1) {
    ldt <- LD[dt$pid,dt$pid]
    diag(ldt) <- 0
    ldt[lower.tri(ldt)] <- 0
    dt$maxld <- apply(abs(ldt),2,max)
    if(thin)
        dt <- dt[maxld<rthr] # ie rsq < 0.01
    dt
}

bdiag_with_dimnames <- function(lst) {
    library(Matrix)
    ret <- Matrix::bdiag(lst)
    rownames(ret) <- colnames(ret) <- unlist(lapply(lst,rownames))
    ret
}
