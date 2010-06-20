###' Create a random-effects module
##'
##' @param bars a list of parsed random-effects terms
##' @param fr a model frame in which to evaluate these terms
##' @param s Number of parameters in the nonlinear mean function (nlmer only)
##'
##' @return a random-effects module that inherits from reTrms
mkReTrms <- function(bars, fr, s = 1L) {
    if (!length(bars))
        stop("No random effects terms specified in formula")
    stopifnot(is.list(bars), all(sapply(bars, is.language)),
              inherits(fr, "data.frame"))
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))

    ## auxiliary {named, for easier inspection}:
    mkBlist <- function(x) {
	ff <- eval(substitute(factor(fac), list(fac = x[[3]])), fr)
	if (all(is.na(ff)))
	    stop("Invalid grouping factor specification, ",
		 deparse(x[[3]]))
	nl <- length(levels(ff))
	mm <- model.matrix(eval(substitute( ~ foo,
					   list(foo = x[[2]]))), fr)
	nc <- ncol(mm)
	nseq <- seq_len(nc)
	sm <- as(ff, "sparseMatrix")
	if (nc	> 1)
	    sm <- do.call(rBind, lapply(nseq, function(i) sm))
	sm@x[] <- t(mm[])
	## When nc > 1 switch the order of the rows of sm
	## so the random effects for the same level of the
	## grouping factor are adjacent.
	if (nc > 1)
	    sm <- sm[as.vector(matrix(seq_len(nc * nl),
				      nc = nl, byrow = TRUE)),]
	list(ff = ff, sm = sm, nl = nl, cnms = colnames(mm))
    }
    blist <- lapply(bars, mkBlist)
    nl <- sapply(blist, "[[", "nl")     # no. of levels per term

    ## order terms stably by decreasing number of levels in the factor
    if (any(diff(nl)) > 0) {
        ord <- rev(order(nl))
        blist <- blist[ord]
        nl <- nl[ord]
    }
    Zt <- do.call(rBind, lapply(blist, "[[", "sm"))
    q <- nrow(Zt)
    ll <- list(Zt = Zt, u = numeric(q))

    ## Create and install Lambda, Lind, etc.  This must be done after
    ## any potential reordering of the terms.
    cnms <- lapply(blist, "[[", "cnms")
    nc <- sapply(cnms, length)          # no. of columns per term
    ncTrms <- nc
    nth <- as.integer((nc * (nc+1))/2)  # no. of parameters per term
    nb <- nc * nl                       # no. of random effects per term
    stopifnot(sum(nb) == q)

    ll$theta <- numeric(sum(nth))

    boff <- cumsum(c(0L, nb))           # offsets into b
    thoff <- cumsum(c(0L, nth))         # offsets into theta
    Lambda <-
        do.call(sparseMatrix,
                do.call(rBind,
                        lapply(seq_along(blist), function(i)
                           {
                               mm <- matrix(seq_len(nb[i]), nc = nc[i],
                                            byrow = TRUE)
                               dd <- diag(nc[i])
                               ltri <- lower.tri(dd, diag = TRUE)
                               ii <- row(dd)[ltri]
                               jj <- col(dd)[ltri]
                               dd[cbind(ii, jj)] <- seq_along(ii)
                               data.frame(i = as.vector(mm[, ii]) + boff[i],
                                          j = as.vector(mm[, jj]) + boff[i],
                                          x = as.double(rep.int(seq_along(ii),
                                          rep.int(nl[i], length(ii))) +
                                          thoff[i]))
                           })))
### FIXME: change this to t(Lambda) before extracting the Lind
    ll$Lind <- as.integer(Lambda@x)

    ## lower bounds on theta elements are 0 if on diagonal, else -Inf
    ll$lower <- -Inf * (ll$theta + 1)
    ll$lower[unique(diag(Lambda))] <- 0
    ll$theta[] <- is.finite(ll$lower)   # initial values of theta are 0 off-diagonal, 1 on
    Lambda@x[] <- ll$theta[ll$Lind]     # initialize elements of Lambda
    ll$Lambda <- Lambda

    ll$Ut <- crossprod(Lambda, Zt)
    ll$Class <- "reTrms"
    if (s > 1) { ## Ut is the sum of vertical sections of Zt
        N <- ncol(Zt)
        ll$Ut <- Reduce("+",
                        lapply(split(seq_len(N),
                                     rep.int(seq_len(s),
                                             rep.int(N %/% s, s))),
                               function(cols) Zt[, cols]))
    }

    ll$L <- Cholesky(tcrossprod(ll$Ut), LDL = FALSE, Imult = 1)
                                        # massage the factor list
    fl <- lapply(blist, "[[", "ff")
    asgn <- seq_along(fl)
                                        # check for repeated factors
    fnms <- names(fl)
    if (length(fnms) > length(ufn <- unique(fnms))) {
        fl <- fl[match(ufn, fnms)]
        asgn <- match(fnms, ufn)
    }
    names(fl) <- ufn
    fl <- do.call(data.frame, c(fl, check.names = FALSE))
    attr(fl, "assign") <- asgn
    ll$flist <- fl
    ll$cnms <- cnms
    do.call("new", ll)
} ## {mkReTrms}

##' Create an feModule
##'
##' @param form Formula for the fixed-effects terms
##' @param fr a model frame in which to evaluate the formula
##' @param contrasts a list of constrasts for factors in fr
##' @param reMod the reModule for the model
##' @param sparseX Logical indicator of sparse X
## FIXME: Why not use s = 0 for indication of non-reweightable?
##' @param s Number of columns in the sqrtXwt matrix
##'
##' @return an object that inherits from feModule
mkFeModule <-
    function(form, fr, contrasts, reMod, sparseX, s = 1L) {
                                        # fixed-effects model matrix X
    nb <- nobars(form[[3]])
    if (is.null(nb)) nb <- 1
    form[[3]] <- nb
    X <- model.Matrix(form, fr, contrasts, sparse = sparseX)
    N <- nrow(X)
    p <- ncol(X)
    rownames(X) <- NULL
    ZtX <- reMod@Zt %*% X
    ll <- list(Class = ifelse(sparseX, "spFeMod", "deFeMod"),
               RZX = ZtX,
               X = X,
               beta = numeric(p))
    ll$RX <-
	if (sparseX)
	    ## watch for the case that crossprod(ll$RZX) is more dense than X.X
	    Cholesky(crossprod(X) + crossprod(ll$RZX), LDL = FALSE)
	else
	    chol(crossprod(X))
    do.call("new", ll)
}

mkRespMod <- function(fr, reMod, feMod, family = NULL,
                      nlenv = NULL, nlmod = NULL, checknl = NULL)
{
    n <- nrow(fr)
    N <- nrow(feMod@X)
                                        # components of the model frame
    y <- model.response(fr)
    # avoid problems with 1D arrays, but keep names
    if(length(dim(y)) == 1) {
        nm <- rownames(y)
        dim(y) <- NULL
        if(!is.null(nm)) names(y) <- nm
    }
    weights <- model.weights(fr)
    if (is.null(weights)) weights <- rep.int(1, n)
    else if (any(weights < 0))
        stop(gettext("negative weights not allowed", domain = "R-lme4"))
    offset <- model.offset(fr)
    if (is.null(offset)) offset <- numeric(N)
    if (length(offset) == 1) offset <- rep.int(offset, N)
    else if (length(offset) != N)
        stop(gettextf("number of offsets (%d) should be %d (s * n)",
                      length(offset), N), domain = "R-lme4")
    p <- ncol(feMod@X)
    q <- nrow(reMod@Zt)
    ll <- list(weights = unname(weights), offset = unname(offset))
    if (!is.null(family)) {
        ll$y <- y                       # may get overwritten later
        rho <- new.env()
        rho$etastart <- model.extract(fr, "etastart")
        rho$mustart <- model.extract(fr, "mustart")
        rho$nobs <- n
        if (is.character(family))
            family <- get(family, mode = "function", envir = parent.frame(3))
        if (is.function(family)) family <- family()
        eval(family$initialize, rho)
        family$initialize <- NULL       # remove clutter from str output
        if(is.null(checknl))
            checknl <- !(family$family %in% c("binomial", "poisson"))
        ll$mu <- unname(rho$mustart)
        lr <- as.list(rho)
        ll[names(lr)] <- lr             # may overwrite y, weights, etc.
        ll$weights <- unname(ll$weights)
        ll$y <- unname(ll$y)
        ll$eta <- family$linkfun(ll$mu)
        ll$sqrtrwt <- sqrt(ll$weights/family$variance(ll$mu))
        ll$wtres <- ll$sqrtrwt * (ll$y - ll$mu)
        ll$sqrtXwt <- matrix(ll$sqrtrwt * family$mu.eta(ll$eta))
        ll$family <- family
        ll <- ll[intersect(names(ll), slotNames("glmerResp"))]
        ll$n <- unname(rho$n)           # for the family$aic function
        ll$Class <- "glmerResp"
    } else {
        ll$sqrtrwt <- sqrt(ll$weights)
        ll$y <- unname(as.numeric(y))
        ll$wtres <- numeric(n)
        ll$mu <- numeric(n)
        if (is.null(nlenv)) {
            ll$Class <- "lmerResp"
            ll$REML <- p
            ll$sqrtXwt <- matrix(ll$sqrtrwt)
        } else {
            ll$Class <- "nlmerResp"
            ll$nlenv <- nlenv
            ll$nlmod <- Quote(nlmod)
            eta <- eval(nlmod, nlenv)
            ll$sqrtXwt <- attr(eta, "gradient")
            if (is.null(ll$sqrtXwt))
                stop("The nonlinear model in nlmer must return a gradient attribute")
            ll$pnames <- colnames(ll$sqrtXwt)
        }
        if(is.null(checknl)) checknl <- TRUE # non-glmer case
    }
    if(checknl && is(reMod, "reTrms") &&
       any(n <= (nlev <- sapply(reMod@flist, function(fac) length(levels(fac))))))
	stop("Number of levels of a grouping factor for the random effects\n",
	     "must be less than the number of observations")
    do.call("new", ll)
}

setGeneric("updateDcmp", function(x, dcmp) standardGeneric("updateDcmp"), valueClass = "list")
setMethod("updateDcmp", signature(x = "reTrms", dcmp = "list"),
          function(x, dcmp) {
              dcmp$dims["reTrms"] <- 1L
              dcmp$dims["q"] <- length(x@u)
              dcmp$dims["nth"] <- length(x@theta)
              dcmp
          })
setMethod("updateDcmp", signature(x = "deFeMod", dcmp = "list"),
          function(x, dcmp) {
              dcmp$dims["spFe"] <- 0L
              dcmp$dims["p"] <- length(x@beta)
              dcmp
          })
setMethod("updateDcmp", signature(x = "spFeMod", dcmp = "list"),
          function(x, dcmp) {
              dcmp$dims["spFe"] <- 1L
              dcmp$dims["p"] <- length(x@beta)
              dcmp
          })

.respBase <- function(x, dcmp) {
    n <- length(x@y)
    N <- length(x@offset)
    dcmp$dims["n"]   <- n
    dcmp$dims["N"]   <- N
    dcmp$dims["nmp"] <- n - dcmp$dims["p"]
    dcmp$dims[c("REML", "GLMM", "NLMM")] <- 0L
    dcmp$dims[c("useSc")] <- 1L
    dcmp
}
    
setMethod("updateDcmp", signature(x = "lmerResp", dcmp = "list"),
          function(x, dcmp) {
              dcmp <- .respBase(x, dcmp)
              dcmp$dims["REML"] <- x@REML
              dcmp
          })
setMethod("updateDcmp", signature(x = "glmerResp", dcmp = "list"),
          function(x, dcmp) {
              dcmp <- .respBase(x, dcmp)
              dcmp$dims["GLMM"] <- 1L
              dcmp$dims["useSc"] <- 0L
              dcmp
          })
setMethod("updateDcmp", signature(x = "nlmerResp", dcmp = "list"),
          function(x, dcmp) {
              dcmp <- .respBase(x, dcmp)
              dcmp$dims["NLMM"] <- 1L
              dcmp
          })
setMethod("updateDcmp", signature(x = "nglmerResp", dcmp = "list"),
          function(x, dcmp) {
              dcmp <- .respBase(x, dcmp)
              dcmp$dims["GLMM"] <- 1L
              dcmp$dims["NLMM"] <- 1L
              dcmp$dims["useSc"] <- 0L
              dcmp
          })

.dcmp <- function() {
    nc <- c("ldL2", "ldRX2", "wrss", "ussq", "pwrss", "drsum", "dev",
            "REML", "sigmaML", "sigmaREML")
    nd <- c("N", "n", "nmp", "nth", "p", "q", "nAGQ", "useSc",
            "reTrms", "spFe", "REML", "GLMM", "NLMM")
    list(cmp = structure(rep(NA, length(nc)), .Names = nc),
         dims = structure(rep(NA_integer_, length(nd)), .Names = nd))
}

setMethod("updateDcmp", signature(x = "merMod", dcmp = "list"),
          function (x, dcmp)
              updateDcmp(x@resp, updateDcmp(x@fe, updateDcmp(x@re, dcmp)))
          )
          
S4toEnv <- function(from) {
    stopifnot(isS4(from))
    ## and we want each to assign each of the slots of the slots
    ans <- new.env()
    for (nm in slotNames(from)) {
	sNms <- slotNames(sl <- slot(from, nm))
	if(length(sNms)) # assign each of the slots of sl
	    for (nnm in sNms) assign(nnm, slot(sl, nnm), envir = ans)
	else ## assign sl itself
	    assign(nm, sl, envir = ans)
    }
    ans
}

## Really need a separate  dense / sparse version of this ????
.lmerM2env <- function (from, envclass, compDev = TRUE)
{
    rho <- S4toEnv(from)
    if (exists("REML", envir = rho, inherits = FALSE))
        rho$REML <- rho$REML > 0
    n <- length(rho$y)
    p <- length(rho$beta)
    nc <- sapply(rho$cnms, length)      # no. of columns per term
    if (all(rho$offset == 0)) rho$offset <- numeric(0)
    if (all(rho$weights == 1)) rho$weights <- numeric(0)
    rho$nobs <- n
    rho$p <- p
    rho$nmp <- n - p
    rho$ncTrms <- nc
    rho$diagonalLambda <- all(nc == 1)
    rho$compDev <- compDev
    ## FIXME ? 'verbose' should also be part of one of the class modules (needed in C)
    cl <- rho$call
    ## This is only a cheap guess, possibly incorrect e.g., in 'lmer(*, verbose=verbose)':
    rho$verbose <-
        as.integer(length(v <- as.list(cl)[names(cl) == "verbose"]) &&
                   !identical(FALSE, eval(v)))
    rho$sparseX <- is(from@fe@X, "sparseMatrix")

    ## temporary "workarounds":
    rho$pwrss <- rho$wrss + sum(rho$u ^ 2)
    rho$gamma <- rho$mu ## needed ?

    ## not yet stored (?)
    ## y <- rho$y
    rho$sqrtrWt <- rho$sqrtrwt
    if (all(rho$sqrtrWt == 1)) rho$sqrtrWt <- numeric(0)
    rm(sqrtrwt, envir = rho)

    rho$sqrtXWt <- rho$sqrtXwt
    if (all(rho$sqrtXWt == 1)) rho$sqrtXWt <- numeric(0)
    rm(sqrtXwt, envir = rho)

    ## Fix up some inconsistencies
    rho$Ut <- crossprod(rho$Lambda, rho$Ut)
    rownames(rho$RZX) <- NULL
    asgn <- attr(rho$flist, "assign")
    rho$flist <- do.call(data.frame, rho$flist)
    attr(rho$flist, "assign") <- asgn
    
    sP <- .setPars # see ./lmer.R
    gP <- function() theta
    gB <- function() cbind(lower = lower,
                           upper = rep.int(Inf, length(theta)))
    environment(sP) <- environment(gP) <- environment(gB) <- rho
    new(envclass, setPars = sP, getPars = gP, getBounds = gB)
}
##setAs("lmerMod", "optenv",  function(from) .lmerM2env(from, "optenv"))
##setAs("lmerMod", "lmerenv", function(from) .lmerM2env(from, "lmerenv"))
##setAs("glmerMod", "glmerenv", function(from) .lmerM2env(from, "glmerenv"))
##setAs("glmerMod",   "merenv", function(from) .lmerM2env(from, "glmerenv"))

lmer2 <- function(formula, data, REML = TRUE, sparseX = FALSE,
                  control = list(), start = NULL,
                  verbose = 0, doFit = TRUE, compDev = TRUE,
                  ## TODO: optimizer = c("bobyqa", "nlminb", "optimize", "optim"),
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
    dcmp <- updateDcmp(reTrms, .dcmp())
                                        # fixed-effects module
    feMod <- mkFeModule(formula, fr, contrasts, reTrms, sparseX)
    respMod <- mkRespMod(fr, reTrms, feMod)
    if (!REML) respMod@REML <- 0L
    ans <- new("merMod",
               call = mc,
               devcomp = updateDcmp(respMod, updateDcmp(feMod, dcmp)),
               frame = fr,
	       re = reTrms, fe = feMod, resp = respMod)
    if (doFit) {                        # optimize estimates
        if (verbose) control$iprint <- 2L
        devfun <- mkdevfun(ans, compDev = compDev)
        opt <- bobyqa(ans@re@theta, devfun, ans@re@lower, control = control)
        ans <- updateMod(ans, opt$par, opt$fval)
    }
    ans
}

updateMod <- function(mod, pars, fval) {
    stopifnot(is(mod, "merMod"),
              is.numeric(pars),
              is.numeric(fval),
              is.numeric(u <- attr(fval, "u")),
              is.numeric(beta <- attr(fval, "beta")))
    re <- mod@re
    fe <- mod@fe
    resp <- mod@resp
    devcomp <- mod@devcomp
    lower <- re@lower
    stopifnot(length(pars) >= length(lower),
              length(beta) == ncol(fe@X),
              length(u) == nrow(re@Zt))
    pars <- pars[seq_along(lower)]
    lst <- .Call(updateDc, mod, pars, beta, u)
    lst[ifelse(is(resp, "lmerResp") && resp@REML, "REML", "dev")] <- as.vector(fval)
    lst["wrss"] <- wrss <- sum(lst[["wtres"]]^2)
    lst["ussq"] <- ussq <- sum(u^2)
    lst["pwrss"] <- pwrss <- wrss + ussq
    lst["sigmaML"] <- sqrt(pwrss/length(resp@y))
    lst["sigmaREML"] <- sqrt(pwrss/(length(resp@y)-length(beta)))
    dc <- mod@devcomp
    nms <- intersect(names(lst), names(dc$cmp))
    dc$cmp[nms] <- unlist(lst[nms])
    mod@re@theta <- pars
    mod@fe@beta <- beta
    mod@re@u <- u
    mod@re@L <- lst[["L"]]
    mod@re@Lambda <- lst[["Lambda"]]
    mod@fe@RX <- lst[["RX"]]
    mod@fe@RZX <- lst[["RZX"]]
    mod@resp@mu <- lst[["mu"]]
    mod@devcomp <- dc
    mod
}
              
##' <description>
##' Returns a function that evaluates the deviance from parameter values
##' <details>
##' From an merMod object create an R function that takes a single
##' argument, which is the new parameter value, and returns the
##' deviance
##' @title Create a deviance evaluation function from an merMod object
##' @param mod an object that inherits from class merMod
##' @param nAGQ number of points per axis for adaptive Gauss-Hermite
##' quadrature. 0 indicates the PIRLSBeta algorithm, 1 is the Laplace approximation
##' @param u0 starting estimate for the PIRLS algorithm
##' @param compDev for lmerMod objects, should the compiled deviance
##' evaluation be used?  Setting compDev to FALSE provides an
##' evaluation using functions from the Matrix package only.
##' @return a function of one argument
##' @author Douglas Bates
mkdevfun <- function(mod, nAGQ = 1L, u0 = numeric(length(mod@re@u)), verbose = 0L, compDev = TRUE) {
    stopifnot(is(mod, "merMod"))
    resp <- mod@resp
    beta0 <- numeric(length(mod@fe@beta))
    if (is(resp, "lmerResp")) {
        if (compDev) return(function(th) .Call(merDeviance, mod, th, beta0, u0, 0L, 3L))
        WtMat <- Diagonal(x = resp@sqrtrwt)
        return(function(th) {
            lower <- mod@re@lower
            stopifnot(is.numeric(th),
                      length(th) == length(lower),
                      all(lower <= th))
            Lambda <- mod@re@Lambda
            Lambda@x <-th[mod@re@Lind]
            Ut <- mod@re@Zt %*% WtMat
            L <- update(mod@re@L, crossprod(Lambda, Ut), mult = 1)
            r0 <- (resp@y - resp@offset) * resp@sqrtrwt
            Utr <- Ut %*% r0
            V <- WtMat %*% mod@fe@X
            Vtr <- crossprod(V, r0)
            UtV <- Ut %*% V
            cu <- solve(L, solve(L, crossprod(Lambda, Utr), sys = "P"),
                        sys = "L")
            RZX <- solve(L, solve(L, crossprod(Lambda, UtV), sys = "P"),
                         sys = "L")
            if (is(V, "sparseMatrix")) {
                RX <- Cholesky(crossprod(V) - crossprod(RZX))
                beta <- solve(RX, Vtr - crossprod(RZX, cu))@x
            } else {
                RX <- chol(crossprod(V) - crossprod(RZX))
                beta <- solve(RX, solve(t(RX), Vtr - crossprod(RZX, cu)))@x
            }
            u <- solve(L, solve(L, cu - RZX %*% beta, sys = "Lt"), sys = "Pt")
            mu <- (resp@offset + crossprod(mod@re@Zt, Lambda %*% u) + mod@fe@X %*% beta)@x
                                        # penalized, weighted residual sum of squares
            pwrss <- sum(c((resp@y - mu)*resp@sqrtrwt, u@x)^2) 
            ldL2 <- 2 * determinant(L)$mod
            if (resp@REML) {
                nmp <- length(mu) - length(beta)
                ldRX2 <- 2 * determinant(RX)$mod
                ldL2 + ldRX2 + nmp * (1 + log(2 * pi * pwrss/nmp))
            } else {
                n <- length(mu)
                ldL2	     + n   * (1 + log(2 * pi * pwrss/n  ))
            }
        })
    }
    if (nAGQ == 0L)
        return(function(pars) .Call(merDeviance, mod, pars, beta0, u0, verbose, 3L))
    thpars <- seq_along(mod@re@theta)
    function(pars) .Call(merDeviance, mod, pars[thpars], pars[-thpars], u0, verbose, 2L)
}
    
## being brave now:
lmer <- lmer2

setMethod("simulate", "merMod",
	  function(object, nsim = 1, seed = NULL, use.u = FALSE, ...)
      {
          stopifnot((nsim <- as.integer(nsim[1])) > 0,
                    is(x, "merMod"), is(x@resp, "lmerResp"))
	  if(!is.null(seed)) set.seed(seed)
	  if(!exists(".Random.seed", envir = .GlobalEnv))
	      runif(1) # initialize the RNG if necessary

	  n <- nrow(X <- object@fe @ X)
	  ## result will be matrix  n x nsim :
	  as.vector(X %*% object@fe @ beta) +  # fixed-effect contribution
	      sigma(object) * (## random-effects contribution:
			       if(use.u) {
				   object@re @ u
			       } else {
				   U <- t(object@re @ Ut)
				   q <- ncol(U)
				   as(U %*% matrix(rnorm(q * nsim), nc = nsim), "matrix")
			       }
			       ## residual contribution:
			       + matrix(rnorm(n * nsim), nc = nsim))
      })

### bootMer() --- <==>  (TODO: semi-*)parametric bootstrap
### -------------------------------------------------------
## Doc: show how  this is equivalent - but faster than
##              boot(*, R = nsim, sim = "parametric", ran.gen = simulate(.,1,.), mle = x)
## --> return a "boot" object -- so we could use boot.ci() etc
## TODO: also allow "semi-parametric" model-based bootstrap:
##    resampling the (centered!) residuals (for use.u=TRUE) or for use.u=FALSE,
##    *both* the centered u's + centered residuals
##    instead of using  rnorm()

##' <description>
##'  Perform model-based (Semi-)parametric bootstrap for mixed models;
##'  The working name for bootMer() was simulestimate(), but this is serious..
##' <details>
##'  ...
##' @title Model-based (Semi-)Parametric Bootstrap for Mixed Models
##' @param x fitted *lmer() model
##' @param FUN a function(x) computing the *statistic* of interest,
##' which must be a numeric vector, possibly named.
##' @param nsim number of simulations, positive integer; the bootstrap
##' 'B' (or 'R').
##' @param seed optional argument to \code{\link{set.seed}}.
##' @param use.u
##' @param verbose logical indicating if progress should print output
##' @param control
##' @return an object of S3 class "boot" {compatible with boot package}
##' @author Martin Maechler
bootMer <- function(x, FUN, nsim = 1, seed = NULL, use.u = FALSE,
                    verbose = FALSE, control = list())
{
    stopifnot((nsim <- as.integer(nsim[1])) > 0,
              is(x, "merMod"), is(x@resp, "lmerResp"))
    FUN <- match.fun(FUN)
    if(!is.null(seed)) set.seed(seed)
    else if(!exists(".Random.seed", envir = .GlobalEnv))
        runif(1) # initialize the RNG if necessary

    mc <- match.call()
    t0 <- FUN(x)
    if (!is.numeric(t0))
        stop("bootMer currently only handles functions that return numeric vectors")

    ## simplistic approach {original for old-lme4 by DB was much smarter}

    n <- nrow(X <- x@fe @ X)
    if(use.u) {
        u <- x@re @ u
    } else {
        U <- t(x@re @ Ut)
        q <- ncol(U)
    }
    Zt <- x@re @ Zt
    X.beta <- as.vector(X %*% x@fe @ beta) # fixed-effect contribution
    sigm.x <- sigma(x)

    ## Here, and below ("optimize"/"bobyqa") using the "logic" of lmer2() itself:
## lmer..Update <- if(is(x, "lmerSp")) lmerSpUpdate else lmerDeUpdate
    devfun <- mkdevfun(x)
##    oneD <- length(x@re@theta) < 2
    theta0 <- x@re@theta
    ## just for the "boot" result -- TODOmaybe drop
    mle <- list(beta = x@fe @ beta, theta = theta0, sigma = sigm.x)

    t.star <- matrix(t0, nr = length(t0), nc = nsim)
    for(i in 1:nsim) {
        y <- {
            X.beta + sigm.x *
                ((if(use.u) u else as.vector(U %*% rnorm(q))) + rnorm(n))
            ##      random effects  contribution            +     Error
        }
	x @ resp @ y <- y

        ## if (oneD) { # use optimize
        ##     d0 <- devfun(0)
        ##     opt <- optimize(devfun, c(0, 10))
        ##     ##                      -------- <<< arbitrary
        ##     ## FIXME ?! if optimal theta > 0, optimize will *not* warn!
        ##     if (d0 <= opt$objective) ## prefer theta == 0 when close
        ##         devfun(0) # -> theta  := 0  and update the rest
        ## } else {
        bobyqa(theta0, devfun, x@re@lower, control = control)
            ## FIXME: also here, prefer \hat\sigma^2 == 0 (exactly)
##        }
        foo <- tryCatch(FUN(x), error = function(e)e)
        if(verbose) { cat(sprintf("%5d :",i)); str(foo) }
        t.star[,i] <- if (inherits(foo, "error")) NA else foo
    }
    rownames(t.star) <- names(t0)

## boot() ends with the equivalent of
    ## structure(list(t0 = t0, t = t.star, R = R, data = data, seed = seed,
    ##		      statistic = statistic, sim = sim, call = call,
    ##		      ran.gen = ran.gen, mle = mle),
    ##		 class = "boot")
    structure(list(t0 = t0, t = t(t.star), R = nsim, data = x@frame,
		   seed = .Random.seed,
		   statistic = FUN, sim = "parametric", call = mc,
		   ## these two are dummies
		   ran.gen = "simulate(<lmerMod>, 1, *)", mle = mle),
	      class = "boot")
}## {bootMer}


PIRLSest <- function(ans, verbose, control, nAGQ) {
    if (verbose) control$iprint <- 2L
                                        # initial optimization of PLSBeta
    u0 <- numeric(length(ans@re@u))
    opt <- bobyqa(ans@re@theta, mkdevfun(ans, 0L, verbose = verbose), ans@re@lower, control = control)
    if (nAGQ > 0L) {
        u0[] <- attr(opt$fval, "u")
        thpars <- seq_along(ans@re@theta)
        bb <- attr(opt$fval, "beta")
        opt <- bobyqa(c(opt$par, bb), mkdevfun(ans, nAGQ, verbose = verbose),
                      lower = c(ans@re@lower, rep.int(-Inf, length(bb))),
                      control = control)
    }
    updateMod(ans, opt$par, opt$fval)
}

glmer2 <- function(formula, data, family = gaussian, sparseX = FALSE,
                   control = list(), start = NULL, verbose = 0L, nAGQ = 1L,
                   doFit = TRUE, 
                   subset, weights, na.action, offset,
                   contrasts = NULL, mustart, etastart, ...)
{
    verbose <- as.integer(verbose)
    mf <- mc <- match.call()
    if (missing(family)) { ## divert using lmer2()
	mc[[1]] <- as.name("lmer")
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
    dcmp <- updateDcmp(reTrms, .dcmp())
    dcmp$dims[c("nAGQ")] <- nAGQ

    feMod <- mkFeModule(formula, fr, contrasts, reTrms, sparseX, TRUE)
    respMod <- mkRespMod(fr, reTrms, feMod, family)
    ans <- new("merMod",
               call = mc,
               devcomp = updateDcmp(respMod, updateDcmp(feMod, dcmp)),
               frame = fr,
	       re = reTrms, fe = feMod, resp = respMod)
    if (!doFit) return(ans)
    ans <- PIRLSest(ans, verbose, control, nAGQ)
    
}## {glmer2}

## being brave now:
glmer <- glmer2



##' Fit a nonlinear mixed-effects model
##'
##' @param formula a nonlinear mixed model formula (see detailed documentation)
##' @param data an optional data frame containing the variables named in
##'    \code{formula}.  By default the variables are taken from the
##'    environment from which \code{nlmer} is called.
##' @param family 
##' @param start starting estimates for the nonlinear model
##'    parameters, as a named numeric vector
##' @param verbose integer scalar passed to nlminb.  If negative then
##'    diagnostic output from the PIRLS (penalized iteratively
##'    reweighted least squares) step is also provided.
##' @param nAGQ number of adaptive Gauss-Hermite quadrature points to use
##' @param doFit logical scalar.  If FALSE the optimization
##'    environment is returned. Otherwise the parameters are estimated
##'    and an object of S4 class "mer" is returned.
##' @param subset further model specifications as in
##'    \code{\link[stats]{lm}}; see there for details.
##' @param weights  further model specifications as in
##'    \code{\link[stats]{lm}}; see there for details.
##' @param na.action  further model specifications as in
##'    \code{\link[stats]{lm}}; see there for details.
##' @param mustart 
##' @param etastart 
##' @param sparseX 
##' @param contrasts  further model specifications as in
##'    \code{\link[stats]{lm}}; see there for details.
##' @param control a list of control parameters passed to bobyqa.
##' @param ... 

##' @return an object of S4 class "merMod"
nlmer2 <- function(formula, data, family = gaussian, start = NULL,
                   verbose = 0L, nAGQ = 1L, doFit = TRUE,
                   subset,
                   weights, na.action, mustart, etastart, sparseX = FALSE,
                   contrasts = NULL, control = list(), ...)
{
    if (!missing(family)) stop("code not yet written")
    mf <- mc <- match.call()
    m <- match(c("data", "subset", "weights", "na.action",
                 "offset", "etastart", "mustart"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
                                        # check the formula
    formula <- as.formula(formula)
    if (length(formula) < 3) stop("formula must be a 3-part formula")
    nlform <- as.formula(formula[[2]])
                                        # Really do need to check twice
    if (length(nlform) < 3) stop("formula must be a 3-part formula")
    nlmod <- as.call(nlform[[3]])
                                        # check for parameter names in start
    if (is.numeric(start)) start <- list(nlpars = start)
    stopifnot((s <- length(pnames <- names(start$nlpars))) > 0,
              is.numeric(start$nlpars))
    if (!all(pnames %in% (anms <- all.vars(nlmod))))
        stop("not all parameter names are used in the nonlinear model expression")
    fr.form <- nlform
## FIXME: This should be changed to use subbars and subnms.
    fr.form[[3]] <-
        parse(text = paste(setdiff(all.vars(formula), pnames),
                         collapse = ' + '))[[1]]
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())

    ## First create nlenv.  For this the nlpar columns are numeric
    for (nm in pnames) fr[[nm]] <- start$nlpars[[nm]]
    nlenv <- new.env()  # inherit from this environment (or environment(formula)?)
    lapply(all.vars(nlmod),
           function(nm) assign(nm, fr[[nm]], envir = nlenv))

    ## Second, extend the frame and convert the nlpar columns to indicators
    n <- nrow(fr)
    frE <- do.call(rbind, lapply(1:s, function(i) fr)) # rbind s copies of the frame
    for (nm in pnames) # convert these variables in fr to indicators
        frE[[nm]] <- as.numeric(rep(nm == pnames, each = n))
                                        # random-effects module
    reTrms <- mkReTrms(findbars(formula[[3]]), frE, s = s)
    dcmp <- updateDcmp(reTrms, .dcmp())
    dcmp$dims["nAGQ"] <- as.integer(nAGQ)[1]

    fe.form <- nlform
    fe.form[[3]] <- formula[[3]]
    feMod <- mkFeModule(fe.form, frE, contrasts, reTrms,
                        sparseX = sparseX, s = s)
                                        # should this check be in mkFeModule?
    p <- length(feMod@beta)
    if ((qrX <- qr(feMod@X))$rank < p)
        stop(gettextf("rank of X = %d < ncol(X) = %d", qrX$rank, p))
    feMod@beta <- qr.coef(qrX, unlist(lapply(pnames, get, envir = nlenv)))
    respMod <- mkRespMod(fr, reTrms, feMod, nlenv = nlenv, nlmod = nlmod)
    respMod@pnames <- pnames
    ans <- new("merMod",
               call = mc,
               devcomp = updateDcmp(respMod, updateDcmp(feMod, dcmp)),
               frame = fr, re = reTrms, fe = feMod, resp = respMod)
    if (!doFit) return(ans)
    PIRLSest(ans, verbose, control, nAGQ)
}

## Methods for the merMod class
setMethod("fixef",  "merMod", function(object, ...)
          structure(object@fe@beta, names = dimnames(object@fe@X)[[2]]))
 
setMethod("formula",  "merMod", function(x, ...) formula(x@call, ...))

##' Extract the random effects.
##'
##' Extract the conditional modes, which for a linear mixed model are
##' also the conditional means, of the random effects, given the
##' observed responses.  These also depend on the model parameters.
##'
##' @param object an object that inherits from the \code{\linkS4class{mer}} class
##' @param postVar logical scalar - should the posterior variance be returned
##' @param drop logical scalar - drop dimensions of single extent
##' @param whichel - vector of names of factors for which to return results

##' @return a named list of arrays or vectors, aligned to the factor list

setMethod("ranef", signature(object = "merMod"),
          function(object, postVar = FALSE, drop = FALSE,
                   whichel = names(ans), ...)
      {
          re <- object@re
          ans <- as.vector(re@Lambda %*% re@u)
          if (is(re, "reTrms")) {
              ## evaluate the list of matrices
              levs <- lapply(fl <- re@flist, levels)
              asgn <- attr(fl, "assign")
              cnms <- re@cnms
              nc <- sapply(cnms, length)
              nb <- nc * (nl <- unlist(lapply(levs, length))[asgn])
              nbseq <- rep.int(seq_along(nb), nb)
              ml <- split(ans, nbseq)
              for (i in seq_along(ml))
                  ml[[i]] <- matrix(ml[[i]], nc = nc[i], byrow = TRUE,
                                    dimnames = list(NULL, cnms[[i]]))
              ## create a list of data frames corresponding to factors
              ans <- lapply(seq_along(fl),
                            function(i)
                            data.frame(do.call(cbind, ml[asgn == i]),
                                       row.names = levs[[i]],
                                       check.names = FALSE))
              names(ans) <- names(fl)
                                        # process whichel
              stopifnot(is(whichel, "character"))
              whchL <- names(ans) %in% whichel
              ans <- ans[whchL]
              
              if (postVar) {
                  vv <- .Call(reTrmsCondVar, re, sigma(object))
                  for (i in seq_along(ans))
                      attr(ans[[i]], "postVar") <- vv[[i]]
              }
              if (drop)
                  ans <- lapply(ans, function(el)
                            {
                                if (ncol(el) > 1) return(el)
                                pv <- drop(attr(el, "postVar"))
                                el <- drop(as.matrix(el))
                                if (!is.null(pv))
                                    attr(el, "postVar") <- pv
                                el
                            })
              class(ans) <- "ranef.mer"
          }
          ans
      })

setMethod("sigma", "merMod", function(object, ...)
      {
          dc <- object@devcomp
          dd <- dc$dims
          if (!dd["useSc"]) return(1.)
          unname(dc$cmp[ifelse(dd["REML"], "sigmaREML", "sigmaML")])
      })

setMethod("model.matrix", signature(object = "merMod"),
	  function(object, ...) object@fe@X)

setMethod("terms", signature(x = "merMod"),
	  function(x, ...) attr(x@frame, "terms"))

setMethod("model.frame", signature(formula = "merMod"),
	  function(formula, ...) formula@frame)

setMethod("deviance", signature(object="merMod"),
	  function(object, REML = NULL, ...) {
              if (!missing(REML)) stop("REML argument not supported")
              unname(object@devcomp$cmp["dev"])
          })

setMethod("logLik", signature(object="merMod"),
	  function(object, REML = NULL, ...)
      {
          if (!missing(REML)) stop("REML argument not supported")
          dc <- object@devcomp
          dims <- dc$dims
          val <- - unname(dc$cmp["dev"])/2
          attr(val, "nall") <- attr(val, "nobs") <- unname(dims["n"])
          attr(val, "df") <- unname(dims["p"] + dims["nth"] + dims["useSc"])
          class(val) <- "logLik"
          val
      })

setMethod("update", signature(object = "merMod"), updateMer)

setMethod("print", "merMod", printMerenv)
setMethod("show",  "merMod", function(object) printMerenv(object))

setMethod("coef", signature(object = "merMod"), coefMer)

setMethod("devcomp", "merMod", function(x, ...) x@devcomp)

setMethod("getL", "reModule", function(x) x@L)
setMethod("getL", "merMod", function(x) x@re@L)

setMethod("isREML", "merMod", function(x) as.logical(x@devcomp$dims["REML"]))

setMethod("refitML", "merMod",
          function (x) {
              if (!isREML(x)) return(x)
              x@fe@REML <- 0L
              x@devcmp$dims["REML"] <- 0L
              bobyqa(x@re@theta, mkdevfun(x), x@re@lower)
              updateMod(x)
          })

setMethod("getCall", "merMod",	function(x) x@call)

setMethod("anova", signature(object = "merMod"), anovaLmer)

setMethod("vcov", signature(object = "merMod"),
	  function(object, correlation = TRUE, sigm = sigma(object), ...)
	  mkVcov(sigm, RX = object@fe@RX, nmsX = colnames(object@fe@X),
		 correlation=correlation, ...))

setMethod("VarCorr", signature(x = "merMod"),
          function(x)
      {
          re <- x@re
          if (!is(re, "reTrms"))
              stop("VarCorr methods require reTrms, not just reModule")
          cnms <- re@cnms
          nc <- sapply(cnms, length)      # no. of columns per term
          mkVarCorr(sigma(x), cnms=cnms, nc = nc,
                    theta=re@theta, flist=re@flist)
      })

setMethod("summary", "merMod",
          function(object, varcov = FALSE, ...)
      {
          resp <- object@resp
          devC <- object@devcomp
          dd <-devC$dims
          cmp <- devC$cmp
          ## FIXME: You can't count on re@flist unless is(re, "reTrms")
          flist <- object@re@flist
          useSc <- as.logical(dd["useSc"])
          sig <- sigma(object)
          REML <- isREML(object)

          fam <- NULL
          if(is(resp, "glmerResp")) fam <- resp@family
          coefs <- cbind("Estimate" = fixef(object),
                         "Std. Error" = sig * sqrt(unscaledVar(RX = object@fe@RX)))
          if (nrow(coefs) > 0) {
              coefs <- cbind(coefs, coefs[,1]/coefs[,2], deparse.level=0)
              colnames(coefs)[3] <- paste(if(useSc) "t" else "z", "value")
          }
          mName <- paste(switch(1L + dd["GLMM"] * 2L + dd["NLMM"],
                                "Linear", "Nonlinear",
                                "Generalized linear", "Generalized nonlinear"),
                         "mixed model fit by",
                          ifelse(REML, "REML", "maximum likelihood"))
          llik <- logLik(object)   # returns NA for a REML fit - maybe change?
          AICstats <- {
              if (REML) cmp["REML"] # do *not* show likelihood stats here
              else {
                  c(AIC = AIC(llik), BIC = BIC(llik), logLik = c(llik),
                    deviance = deviance(object))
              }
          }
          varcor <- VarCorr(object)
                                        # use S3 class for now
          structure(list(methTitle = mName,
                         devcomp = devC, isLmer = is(resp, "lmerResp"), useScale = useSc,
                         logLik = llik, family = fam,
                         ngrps = sapply(flist, function(x) length(levels(x))),
                         coefficients = coefs,
                         sigma = sig,
                         vcov = if(varcov) vcov(object),
                         varcor = varcor, # and use formatVC(.) for printing.
                         AICtab= AICstats,
                         call = object@call
                         ), class = "summary.merenv")
      })
