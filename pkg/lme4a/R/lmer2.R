##' Create a random-effects module
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
    ll <- list(Zt = Zt, u = numeric(q), Utr = numeric(q))

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
    ll$ldL2 <- numeric(1)
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
    ll$nlev <- nl
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
               UtV = ZtX,
               V = X,
               VtV = crossprod(X),
               Vtr = numeric(p),
               X = X,
               beta = numeric(p),
               ldRX2= numeric(1))
    ll$RX <-
	if (sparseX)
	    ## watch for the case that crossprod(ll$RZX) is more dense than X.X
	    Cholesky(ll$VtV + crossprod(ll$RZX), LDL = FALSE)
	else
	    chol(ll$VtV)
    if (s > 1) {
        ll$V <- Reduce("+", lapply(split(seq_len(N) ,
                                         rep.int(1:s, rep.int(N %/% s, s))),
                                   function(rows) X[rows, ]))
    }
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
    ll <- list(weights = unname(weights), offset = unname(offset), wrss = numeric(1))
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
        ll$muEta <- family$mu.eta(ll$eta)
        ll$var <- family$variance(ll$mu)
        ll$sqrtrwt <- sqrt(ll$weights/ll$var)
        ll$wtres <- ll$sqrtrwt * (ll$y - ll$mu)
        ll$sqrtXwt <- matrix(ll$sqrtrwt * ll$muEta)
        ll$family <- family
        ll$devres <- numeric(1)
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
    if(checknl && any(reMod@nlev >= n))
	stop("Number of levels of a grouping factor for the random effects\n",
	     "must be less than the number of observations")
    do.call("new", ll)
}


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
    n <- length(rho$y)
    p <- length(rho$beta)
    nc <- sapply(rho$cnms, length)      # no. of columns per term
    rho$nobs <- n
    rho$nmp <- n - p
    rho$ncTrms <- nc
    rho$diagonalLambda <- all(nc == 1)
    rho$compDev <- compDev
    ## FIXME ? 'verbose' should also be part of one of the class modules (needed in C)
    cl <- rho$call
    ## This is only a cheap guess, possibly incorrect e.g., in 'lmer(*, verbose=verbose)':
    rho$verbose <-
        length(v <- as.list(cl)[names(cl) == "verbose"]) && !identical(FALSE, eval(v))
    rho$sparseX <- is(from@fe@X, "sparseMatrix")

    ## temporary "workarounds":
    rho$pwrss <- rho$wrss + sum(rho$u ^ 2)
    rho$gamma <- rho$mu ## needed ?

    ## not yet stored (?)
    ## y <- rho$y
    rho$sqrtrWt <- sqrt(rho$weights)
    ##  sqrtXWt  :=  sqrt of model matrix row weights
    rho$sqrtXWt <- sqrt(rho$weights)    # to ensure a distinct copy

    sP <- .setPars # see ./lmer.R
    gP <- function() theta
    gB <- function() cbind(lower = lower,
                           upper = rep.int(Inf, length(theta)))
    environment(sP) <- environment(gP) <- environment(gB) <- rho
    new(envclass, setPars = sP, getPars = gP, getBounds = gB)
}
setAs("lmerMod", "optenv",  function(from) .lmerM2env(from, "optenv"))
setAs("lmerMod", "lmerenv", function(from) .lmerM2env(from, "lmerenv"))
setAs("glmerMod", "glmerenv", function(from) .lmerM2env(from, "glmerenv"))
setAs("glmerMod",   "merenv", function(from) .lmerM2env(from, "glmerenv"))

lmer2 <- function(formula, data, REML = TRUE, sparseX = FALSE,
                  control = list(), start = NULL,
                  verbose = 0, doFit = TRUE,
                  ## think of also copying 'compDev = FALSE' from lmer1()
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
                                        # fixed-effects module
    feMod <- mkFeModule(formula, fr, contrasts, reTrms, sparseX)
    respMod <- mkRespMod(fr, reTrms, feMod)
    reTrms@Utr[] <- (reTrms@Zt %*% respMod@y)@x
    feMod@Vtr[] <- crossprod(feMod@X, respMod@y)@x
    if (!REML) respMod@REML <- 0L
    ans <- new(ifelse (sparseX, "lmerSp", "lmerDe"), call = mc, frame = fr,
	       re = reTrms, fe = feMod, resp = respMod)
    if (doFit) {                        # optimize estimates
        code <- if(is(ans, "lmerSp")) lmerSpUpdate else lmerDeUpdate
	if(verbose) {
	    ..it <- 0L
	    f.width <- (f.dig <- getOption("digits")) + 5 # "+5": e.g. for 1.2e-05
	    devfun <- function(th) {
		r <- .Call(code, ans, th)
		..it <<- ..it + 1L
		cat(sprintf("%3d : %s |-> %15.12g\n", ..it,
			    paste(sapply(th, format, width = f.width),
                                  collapse = " "),
                            r))
		r
	    }
	} else
	    devfun <- function(th) .Call(code, ans, th)

        if (length(ans@re@theta) < 2) { # use optimize
            d0 <- devfun(0)
            opt <- optimize(devfun, c(0, 10))
            ##                      -------- <<< arbitrary
            ## FIXME ?! if optimal theta > 10, optimize will *not* warn!
            if (d0 <= opt$objective) { ## prefer theta == 0 when close
                cat(sprintf("dev(th =0) = %.12g <= %.12g = opt$obj(th=%g)%s\n",
                            d0, opt$objective, opt$minimum, " --> th := 0"))
                devfun(0) # -> theta  := 0  and update the rest
            }
        } else {
            ## if (verbose > 1) control$iprint <- verbose
            bobyqa(ans@re@theta, devfun, ans@re@lower, control = control)
            ## FIXME: also here, prefer \hat{\sigma_a^2} == 0 (exactly)
        }
    }
    ans
}

## being brave now:
lmer <- lmer2

setMethod("simulate", "lmerMod",
	  function(object, nsim = 1, seed = NULL, use.u = FALSE, ...)
      {
          stopifnot((nsim <- as.integer(nsim[1])) > 0) ## is(x, "lmer")
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
##' @param FUN a function(x) computating the *statistic* of interest,
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
              is(x, "lmerMod"))
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
    lmer..Update <- if(is(x, "lmerSp")) lmerSpUpdate else lmerDeUpdate
    devfun <- function(th) .Call(lmer..Update, x, th)
    oneD <- length(x@re@theta) < 2

    theta0 <- x@re@theta # to use as starting value
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
	x @ fe @ Vtr <- crossprod(X, y)@x
	x @ re @ Utr <- (Zt %*% y)@x

        if (oneD) { # use optimize
            d0 <- devfun(0)
            opt <- optimize(devfun, c(0, 10))
            ##                      -------- <<< arbitrary
            ## FIXME ?! if optimal theta > 0, optimize will *not* warn!
            if (d0 <= opt$objective) ## prefer theta == 0 when close
                devfun(0) # -> theta  := 0  and update the rest
        } else {
            bobyqa(theta0, devfun, x@re@lower, control = control)
            ## FIXME: also here, prefer \hat\sigma^2 == 0 (exactly)
        }
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

glmer2 <- function(formula, data, family = gaussian, sparseX = FALSE,
                   control = list(), start = NULL, verbose = 0L, doFit = TRUE,
                   subset, weights, na.action, offset,
                   contrasts = NULL, nAGQ = 1, mustart, etastart, ...)
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

    if (!(all(1 == (nAGQ <- as.integer(nAGQ)))))
        warning("nAGQ > 1 has not been implemented, using Laplace")
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
    reTrms <- mkReTrms(findbars(formula[[3]]), fr) # random-effects module
    feMod <- mkFeModule(formula, fr, contrasts, reTrms, sparseX, TRUE)
    respMod <- mkRespMod(fr, reTrms, feMod, family)
    feMod@V <- Diagonal(x = respMod@sqrtXwt[,1]) %*% feMod@X
    ans <- new(ifelse(sparseX, "glmerSp", "glmerDe"), call = mc,
               frame = fr, re = reTrms, fe = feMod, resp = respMod)
    .Call(glmerDeIRLS, ans, verbose)
    if (doFit) {                        # optimize estimates
        code <- if(is(ans, "glmerSp")) glmerSpPIRLSBeta else glmerDePIRLSBeta
	if(verbose) {
	    ..it <- 0L
	    f.width <- (f.dig <- getOption("digits")) + 5 # "+5": e.g. for 1.2e-05
	    devfun <- function(th) {
		r <- .Call(code, ans, th, verbose)
		..it <<- ..it + 1L
		cat(sprintf("%3d : %s |-> %15.12g\n", ..it,
			    paste(sapply(th, format, width = f.width),
                                  collapse = " "),
                            r))
		r
	    }
	} else
	    devfun <- function(th) .Call(code, ans, th, verbose)

        if (length(ans@re@theta) < 2) { # use optimize
            d0 <- devfun(0)
            opt <- optimize(devfun, c(0, 10))
            ##                      -------- <<< arbitrary
            ## FIXME ?! if optimal theta > 10, optimize will *not* warn!
            if (d0 <= opt$objective) { ## prefer theta == 0 when close
                cat(sprintf("dev(th =0) = %.12g <= %.12g = opt$obj(th=%g)%s\n",
                            d0, opt$objective, opt$minimum, " --> th := 0"))
                devfun(0) # -> theta  := 0  and update the rest
            }
        } else {
            ## if (verbose > 1) control$iprint <- verbose
            bobyqa(ans@re@theta, devfun, ans@re@lower, control = control)
            ## FIXME: also here, prefer \hat{\sigma_a^2} == 0 (exactly)
        }

        ## now PIRLS optimization to refine the answer
        code <- if(is(ans, "glmerSp")) glmerSpPIRLS else glmerDePIRLS
        thpars <- seq_along(ans@re@theta)
        bb <- ans@fe@beta
        devfun <- function(pars) {
            .Call(feSetBeta, ans@fe, pars[-thpars])
            .Call(code, ans, pars[thpars], verbose)
        }
        if (verbose) control$iprint <- 2L
        bobyqa(c(ans@re@theta, bb), devfun,
               c(ans@re@lower, rep.int(-Inf, length(bb))),
               control = control)
        .Call(glmerDeUpdateRzxRx, ans)
    }
    ans
}## {glmer2}

## being brave now:
glmer <- glmer2



##' Fit a nonlinear mixed-effects model
##'
##' @param formula a nonlinear mixed model formula (see detailed documentation)
##' @param data an optional data frame containing the variables named in
##'    \code{formula}.  By default the variables are taken from the
##'    environment from which \code{nlmer} is called.
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
##' @param contrasts  further model specifications as in
##'    \code{\link[stats]{lm}}; see there for details.
##' @param control a list of control parameters passed to nlminb.  The
##'    defaults are given in the (hidden) function \code{lmerControl}.

##' @return if doFit is FALSE an environment, otherwise an object of S4 class "mer"
nlmer2 <- function(formula, data, family = gaussian, start = NULL,
                   verbose = 0, nAGQ = 1, doFit = TRUE, subset,
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

    fe.form <- nlform
    fe.form[[3]] <- formula[[3]]
    feMod <- mkFeModule(fe.form, frE, contrasts, reTrms,
                        sparseX = sparseX, s = s)
                                        # should this check be in mkFeModule?
    p <- length(feMod@beta)
    if ((qrX <- qr(feMod@X))$rank < p)
        stop(gettextf("rank of X = %d < ncol(X) = %d", qrX$rank, p))
    feMod@beta[] <- qr.coef(qrX, unlist(lapply(pnames, get, envir = nlenv)))
    respMod <- mkRespMod(fr, reTrms, feMod, nlenv = nlenv, nlmod = nlmod)
    respMod@pnames <- pnames
    ans <- new(ifelse(sparseX, "nlmerSp", "nlmerDe"), call = mc,
               frame = fr, re = reTrms, fe = feMod, resp = respMod)
    ans
}

## #    rho$eta <- numeric(length(eta))
## #    rho$eta[] <- eta

##     rho$start <- numeric(p)  # must be careful that these are distinct
##     rho$start[] <-
##     rho$RX <- qr.R(qrX)

##     q <- length(rho$u)
##     rho$u0 <- numeric(q)
##     if (!missing(verbose)) control$trace <- as.integer(verbose[1])
##     rho$dims["verb"] <- control$trace
## #    control$trace <- abs(control$trace) # negative values give PIRLS output
##     rho$control <- control
##     rho$mc <- mc                        # store the matched call

##     rho$bds <- getBounds(rho)
##     theta0 <- getPars(rho)
##     if (!is.null(start$theta)) {
##         stopifnot(length(st <- as.double(start$theta)) == length(theta0))
##         setPars(rho, st)
##     } else {
##         setPars(rho, theta0) # one evaluation to check structure
##     }
##     rho$beta0[] <- rho$beta
##     rho$u0[] <- rho$u
##     if (!doFit) return(rho)

## #    merFinalize(rho)
## }## {nlmer}
