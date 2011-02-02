##' Solve for both sets of coefficients
##'
##' @title Solve for both beta and u
##' @param rem random-effects module
##' @param fem fixed-effects module
##' @param Xwts square root of the X weights
##' @param resids weighted residuals
##' @param useU0 should the gradient be based on the value of u0 or u
##' @return none
solveBetaU <- function(rem, fem, Xwts, resids) {
### FIXME: change the Xwts to be a vector, not a matrix
    Xwts <- as.matrix(Xwts)
    rem$reweight(Xwts, resids, TRUE)    # update U, Utr and L
    fem$reweight(Xwts, resids)          # update V, VtV and Vtr
    fem$updateUtV(rem$Ut)
    fem$updateRzxpRxpp(rem$Lambdap, rem$Lp)
    tmp <- fem$updateIncr(rem$cu)       # increment for beta
    rem$updateIncr(tmp)                 # increment for u
}

lmer2 <- function(formula, data, REML = TRUE, sparseX = FALSE,
                  control = list(), start = NULL,
                  verbose = 0L, doFit = TRUE, compDev = TRUE,
                  subset, weights, na.action, offset,
                  contrasts = NULL, gctort=FALSE, ...)
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
        solveBetaU(rem, fem, resp$sqrtXwt, resp$wtres)
                                        # update mu (full step)
        resp$updateMu(rem$linPred1(1) + fem$linPred1(1))
                                        # profiled deviance or REML
        resp$Laplace(rem$ldL2, fem$ldRX2, rem$sqrLenU)
    }
    ## one evaluation of devfun to ensure that all values are set
if (gctort) gctorture(1)
    opt <- list(fval=devfun(reTrms@theta))

    if (doFit) {                        # optimize estimates
        if (verbose) control$iprint <- as.integer(verbose)
        opt <- bobyqa(reTrms@theta, devfun, reTrms@lower, control = control)
    }
if (gctort) gctorture(0)
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

##' Create an lmerResp, glmerResp or (later) nlmerResp instance
##' 
##' @title Create a [ng]lmerResp instance
##' @param fr a model frame
##' @param family the optional glm family (glmRespMod only)
##' @param nlenv the nonlinear model evaluation environment (nlsRespMod only)
##' @param nlmod the nonlinear model function (nlsRespMod only)
##' @return a lmerResp (or glmerResp) instance
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
stepFac <- function(rem, fem, resp, verbose=FALSE) {
    pwrss0 <- resp$wrss + rem$sqrLenU
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
}

pwrssUpdate <- function(rem, fem, resp, verbose) {
    for (i in 1:5) {
        solveBetaU(rem, fem, resp$sqrtXwt, resp$wtres)
        stepFac(rem, fem, resp, verbose)
    }
}

pwrssUpdate2 <- function(rem, fem, resp, verbose) {
    repeat {
        rem$reweight(resp$sqrtXwt, resp$wtres, TRUE) 
        ccrit <- sqrt(rem$solveIncr()/resp$pwrss)
        if (ccrit < 0.001) break
        stepFac(rem, fem, resp, verbose)
    }
}

setMethod("show", signature("Rcpp_glmerResp"), function(object)
      {
          with(object,
               print(head(cbind(weights, offset, eta, mu, y, muEta,
                                variance, sqrtrwt, wtres, sqrtXwt=sqrtXwt[,1],
                                sqrtWrkWt, wrkResids, wrkResp))))
      })

setMethod("show", signature("Rcpp_deFeMod"), function(object)
      {
          with(object, print(cbind(coef0, incr, coef, Vtr)))
      })

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
    fr.form <- subbars(formula) # substitute "|" for "+" -
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
                                        # random-effects module
    reTrms <- mkReTrms(findbars(formula[[3]]), fr)
    rem <- new(reModule, reTrms@Zt, reTrms@Lambda,
               reTrms@L, reTrms@Lind, reTrms@lower)
    rem$theta <- reTrms@theta           # initialize structures
                                        # fixed-effects module
    feMod <- mkFeModule(formula, fr, contrasts, reTrms, sparseX)
    n <- nrow(fr)
    X <- feMod@X
    p <- ncol(X)
    q <- nrow(rem$Zt)
    fem <- new(deFeMod, X, n, p, q)
                                        # response module
    resp <- mkRespMod2(fr, family)
                                        # initial step for working residuals
    fem$reweight(as.matrix(resp$sqrtWrkWt), resp$wrkResp)
    fem$solveIncr()
    resp$updateMu(fem$linPred1(1))
    fem$installCoef0()

    resp$updateWts()
    fem$incr <- numeric(p)
    pwrssUpdate2(rem, fem, resp, verbose)
    resp$updateWts()
    
    resp$pwrss <- resp$updateWts() + rem$sqrLenU
    list(rem = rem, fem = fem, resp = resp)
}## {glmer2}
