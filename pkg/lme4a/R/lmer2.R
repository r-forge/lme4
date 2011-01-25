solveBetaU <- function(rem, fem, Xwts, resids, useU0=FALSE) {
    Xwts <- as.matrix(Xwts)
    rem$reweight(Xwts, resids, useU0)
    fem$reweight(Xwts, resids)
    fem$updateUtV(rem$Ut)                   # update V, VtV and Vtr
    fem$updateRzxpRxpp(rem$Lambdap, rem$Lp) # update RZX and RX
    rem$updateIncr(fem$updateIncr(rem$cu))  # increments to u and beta
}

lmer2 <- function(formula, data, REML = TRUE, sparseX = FALSE,
                  control = list(), start = NULL,
                  verbose = 0, doFit = TRUE, compDev = TRUE,
                  subset, weights, na.action, offset,
                  contrasts = NULL, ...)
{
    mf <- mc <- match.call()
    ## '...' handling up front, safe-guarding against typos ("familiy") :
    if(length(l... <- list(...))) {
        if (!is.null(l...$family)) {  # call glmer if family specified
            mc[[1]] <- as.name("glmer")
            return( eval(mc, parent.frame()) )
        }
        ## Check for method argument which is no longer used
        if (!is.null(method <- l...$method)) {
	    msg <- paste("Argument", sQuote("method"), "is deprecated.")
            if (match.arg(method, c("Laplace", "AGQ")) == "Laplace") {
                warning(msg)
                l... <- l...[names(l...) != "method"]
            } else stop(msg)
        }
        if(length(l...))
            warning("extra arguments ", paste(names(l...), sep=", "),
                    " are disregarded")
    }

    stopifnot(length(formula <- as.formula(formula)) == 3)
    if (missing(data)) data <- environment(formula)
                                        # evaluate and install the model frame :
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula) # substituted "|" by "+" -
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
                                        # random effects and terms modules
    reTrms <- mkReTrms(findbars(formula[[3]]), fr)
##-nL if (any(reTrms@nLevs >= ncol(reTrms@Zt)))
    if (any(unlist(lapply(reTrms@flist, nlevels)) >= ncol(reTrms@Zt)))
        stop("number of levels of each grouping factor must be less than number of obs")
    rem <- new(reModule, reTrms@Zt, reTrms@Lambda,
               reTrms@L, reTrms@Lind, reTrms@lower)
                                        # fixed-effects module
    feMod <- mkFeModule(formula, fr, contrasts, reTrms, sparseX)
    n <- nrow(fr)
    X <- feMod@X
    p <- ncol(X)
    q <- nrow(rem$Zt)
    fem <- new(deFeMod, X, n, p, q)
                                        # response module
    resp <- mkRespMod2(fr)
    if (REML) resp$REML <- ncol(fem$X)
    devfun <- function(theta) {
        rem$theta <- theta              # update theta, Lambda
        resp$updateMu(numeric(n))       # zero the current mu
        solveBetaU(rem, fem, resp$sqrtXwt, resp$wtres, TRUE)
                                        # update mu (full step)
        resp$updateMu(rem$linPred1(1) + fem$linPred1(1))
                                        # profiled deviance or REML
        resp$Laplace(rem$ldL2, fem$ldRX2, rem$sqrLenU)
    }
    devfun(reTrms@theta)                # one evaluation ensure all values are set

    if (doFit) {                        # optimize estimates
        if (verbose) control$iprint <- 2L
        opt <- bobyqa(reTrms@theta, devfun, reTrms@lower, control = control)
    }
    
    sqrLenU <- rem$sqrLenU
    wrss <- resp$wrss
    pwrss <- wrss + sqrLenU
    dims <- c(N=n, n=n, nmp=n-p, nth=length(rem$theta), p=p, q=q,
              nAGQ=NA_integer_, useSc=1L, reTrms=length(reTrms@cnms),
              spFe=0L, REML=resp$REML, GLMM=0L, NLMM=0L)
    cmp <- c(ldL2=rem$ldL2, ldRX2=fem$ldRX2, wrss=wrss, ussq=sqrLenU, pwrss=pwrss, drsum=NA,
             dev=if(REML)NA else opt$fval, REML=if(REML)opt$fval else NA,
             sigmaML=sqrt(pwrss/n), sigmaREML=sqrt(pwrss/(n-p)))
    re <- new("reTrms", flist=reTrms@flist, cnms=reTrms@cnms, L=rem$L, Lambda=rem$Lambda,
              Lind=rem$Lind, Zt=rem$Zt, lower=rem$lower, theta=rem$theta, u=rem$u)
    new("merMod", call=mc, frame=fr, devcomp = list(cmp=cmp,dims=dims), re=re,
        fe=as(fem, "deFeMod"), resp=as(resp, "lmerResp"))
}## { lmer2 }

mkRespMod2 <- function(fr, family = NULL, nlenv = NULL, nlmod = NULL) {
    n <- nrow(fr)
    y <- model.response(fr)
    if(length(dim(y)) == 1) {
        ## avoid problems with 1D arrays, but keep names
        nm <- rownames(y)
        dim(y) <- NULL
        if(!is.null(nm)) names(y) <- nm
    }
    yy <- y
    if (is.factor(yy)) yy <- as.numeric(yy != levels(yy)[1])
    if (is.matrix(yy) && ncol(yy) == 2L) yy <- yy[,1]/rowSums(yy)
    ans <- new(lmerResp, 0L, yy)
    if (!is.null(offset <- model.offset(fr))) {
        if (length(offset) == 1L) offset <- rep.int(offset, n)
        stopifnot(length(offset) == n)
	ans$offset <- unname(offset)
    }
    if (!is.null(weights <- model.weights(fr))) {
        stopifnot(length(weights) == n, all(weights >= 0))
	ans$weights <- unname(weights)
    }
    if (is.null(family)) return(ans)
    rho <- new.env()
    rho$etastart <- model.extract(fr, "etastart")
    rho$mustart <- model.extract(fr, "mustart")
    rho$weights <- ans$weights
    rho$nobs <- n
    eval(family$initialize, rho)
    family$initialize <- NULL       # remove clutter from str output
        
    ans <- new(glmerResp, family, rho$y, rho$weights, ans$offset, rho$n)
    ans$updateMu(family$linkfun(unname(rho$mustart)))
    ans
}

##' Determine a step factor that will reduce the pwrss
##'
##' The penalized, weighted residual sum of squares (pwrss) is the sum
##' of the weighted residual sum of squares from the resp module and
##' the squared length of u from the rem module.  Both the rem and the
##' fem contain a base and an increment for the coefficients.
##' @title Determine a step factor
##' @param rem random-effects module
##' @param fem fixed-effects module
##' @param resp response module
##' @return pwrss under the new weights at the new u0 and coef0
##' @author Douglas Bates
stepFac <- function(rem, fem, resp, verbose=FALSE) {
    pwrss0 <- resp$pwrss
### FIXME: Probably store the squared length of u in the respModule
### and return the penalized, weighted residual sum of squares

### FIXME: It would be trivial to create a sqrLenU0 method for
### reModule, which might make some of the logic cleaner.
    for (fac in 2^(-(0:10))) {
        ## the calculation of pwrss1 is done in two steps because I'm
        ## not sure of the evaluation order and sqrLenU must follow
        ## linPred1(fac)
        pwrss1 <- resp$updateMu(rem$linPred1(fac)+fem$linPred1(fac))
        pwrss1 <- pwrss1 + rem$sqrLenU  
        if (verbose) cat(sprintf("pwrss0=%g, diff=%12.5g, fac=%7.4f\n",
                                 pwrss0, pwrss0 - pwrss1, fac))
        if (pwrss1 < pwrss0) {
            resp$updateWts()
            rem$installU0()
            fem$installCoef0()
            resp$pwrss <- resp$updateMu(rem$linPred1(0)+fem$linPred1(0)) + rem$sqrLenU
            break
        }
        if (fac < 0.001)
            stop("step factor reduced below 0.001 without reducing wrss")
    }
    invisible(NULL)
}


pwrssUpdate <- function(theta, pwrss, rem, fem, resp) {
    rem$theta <- theta
    solveBetaU(rem, fem, resp$sqrtXwt, resp$wtres)
    stepFac(rem, fem, resp)
}

pwrssUpdate2 <- function(theta, rem, fem, resp) {
    rem$theta <- theta
    rem$reweight(resp$sqrtXwt, resp$wtres, FALSE)
    rem$solveIncr()
    stepFac(rem, fem, resp)
}

glmer2 <- function(formula, data, family = gaussian, sparseX = FALSE,
                   control = list(), start = NULL, verbose = 0L, nAGQ = 1L,
                   doFit = TRUE, subset, weights, na.action, offset,
                   contrasts = NULL, mustart, etastart, ...)
{
    verbose <- as.integer(verbose)
    mf <- mc <- match.call()
    if (missing(family)) { ## divert using lmer2()
	mc[[1]] <- as.name("lmer2")
	return(eval(mc, parent.frame()))
    }
### '...' handling up front, safe-guarding against typos ("familiy") :
    if(length(l... <- list(...))) {
	## Check for invalid specifications
	if (!is.null(method <- list(...)$method)) {
	    msg <- paste("Argument", sQuote("method"),
			 "is deprecated.\nUse", sQuote("nAGQ"),
			 "to choose AGQ.  PQL is not available.")
	    if (match.arg(method, c("Laplace", "AGQ")) == "Laplace") {
		warning(msg)
		l... <- l...[names(l...) != "method"]
	    } else stop(msg)
	}
	if(length(l...))
	    warning("extra arguments ", paste(names(l...), sep=", "),
		    " are disregarded")
    }
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if(is.function(family)) family <- family()
    if (family$family %in% c("quasibinomial", "quasipoisson", "quasi"))
        stop('"quasi" families cannot be used in glmer')
    nAGQ <- as.integer(nAGQ)[1]
    if (nAGQ > 1) warning("nAGQ > 1 has not been implemented, using Laplace")
    stopifnot(length(formula <- as.formula(formula)) == 3)
    if (missing(data)) data <- environment(formula)
                                        # evaluate and install the model frame :
    m <- match(c("data", "subset", "weights", "na.action", "offset",
                 "mustart", "etastart"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula) # substituted "|" by "+" -
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
                                        # random-effects module
    reTrms <- mkReTrms(findbars(formula[[3]]), fr)
    rem <- new(reModule, reTrms@Zt, reTrms@Lambda,
               reTrms@L, reTrms@Lind, reTrms@lower)
                                        # fixed-effects module
    feMod <- mkFeModule(formula, fr, contrasts, reTrms, sparseX)
    n <- nrow(fr)
    X <- feMod@X
    p <- ncol(X)
    q <- nrow(rem$Zt)
    fem <- new(deFeMod, X, n, p, q)
                                        # response module
    resp <- mkRespMod2(fr, family)
                                       
    rem$theta <- reTrms@theta           # initialize structures
    solveBetaU(rem, fem, resp$sqrtWrkWt, resp$wrkResp, TRUE)
    resp$updateMu(rem$linPred1(1) + fem$linPred1(1))
    rem$installU0()
    fem$installCoef0()
    resp$pwrss <- resp$updateWts() + rem$sqrLenU
    list(rem = rem, fem = fem, resp = resp)
}## {glmer2}
