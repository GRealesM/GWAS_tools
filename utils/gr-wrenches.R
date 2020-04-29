## Set of useful functions for file processing and other stuff


##' Convert effect sizes and standard error from linear to log-odds scale
##'
##' Binary traits in case-control studies can be analysed by BOLT-LMM, which uses a linear model.
##' This function provides an approximation from linear to log-odds scale for Beta and SE, 
##' as described in BOLT-LMM manual (https://data.broadinstitute.org/alkesgroup/BOLT-LMM/).
##' @title Transform Beta and SE from linear to logg-odds scale
##' @param beta a numeric vector of Beta in linear scale.
##' @param se a numeric vector of SE in linear scale.
##' @param N0 number of controls in the study.
##' @param N1 number of cases in the study.
##' @return a list containing beta and se.beta in log-odds scale.
linORscale <- function(beta, se,N0,N1){
	cp  <- N1/(N0+N1) 
	BETA  <- BETA/(cp * (1-cp))
	SE  <- SE/(cp * (1-cp))
	ORscaled  <- list(BETA,SE)
	return(ORscaled)
}

##' Estimate trait standard deviation given vectors of variance of coefficients, MAF and sample size, and correct BETA and SE.
##' 
##' For studies of quantitative traits, measurements can come in different scales, leading to 
##' differences in standard deviation of the effects.
##' This function estimates the trait standard deviation from MAF and sample size, and by dividing beta and se 
##' by this estimate we can make var.beta = 1, thus making different datasets of quantitative traits in different
##' scales comparable.
##' 
##' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
##' var(X) = 2*maf*(1-maf)
##' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
##' 
##' This is an adaptation of sdY.est function, by Chris Wallace
##'
##' @title Estimate trait variance, internal function
##' @param beta vector of coefficients
##' @param se vector of standard error of coefficients
##' @param maf vector of MAF (same length as beta and se)
##' @param n sample size
##' @return a list containing adjusted Beta and SE by estimated standard deviation of Y
##' @author Guillermo Reales, Chris Wallace
sdY.correction <- function(beta, se, maf, n) {
    if(length(beta) != length(maf) | length(se) != length(maf)) stop("Beta, SE and/or MAF are not the same length")
    oneover <- 1/se^2
    nvx <- 2 * n * maf * (1-maf)
    m <- lm(nvx ~ oneover - 1)
    cf <- coef(m)[['oneover']]
    if(cf < 0)
        stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
    message("Estimated sdY is ", cf)
    BETA  <- beta/sqrt(cf)
    SE  <- se/sqrt(cf)
    adj.BETASE  <- list(BETA,SE)
    return(adj.BETASE)
}


