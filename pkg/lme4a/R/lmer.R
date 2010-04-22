##' Check and install various objects, including frame, y, weights
##' from the model frame in the environment rho.
##'
##' @param fr a model frame
##' @param rho an environment
##' @param muetastart logical - check for mustart and etastart
check_y_weights <- function(fr, rho, muetastart = FALSE) {
    rho$frame <- fr
    attr(rho$frame, "terms") <- NULL
    rho$nobs <- n <- nrow(fr)
                                        # components of the model frame
    y <- model.response(fr)
    # avoid problems with 1D arrays, but keep names
    if(length(dim(y)) == 1) {
        nm <- rownames(y)
        dim(y) <- NULL
        if(!is.null(nm)) names(y) <- nm
    }
    rho$y <- y
    rho$weights <- as.numeric(model.weights(fr))
    if (any(rho$weights < 0))
        stop(gettext("negative weights not allowed", domain = "R-lme4"))
    loff <- length(offset <- as.numeric(model.offset(fr)))
    if (loff) {
        if (loff == 1) {
            offset <- rep.int(offset, n)
        } else if (loff != n) {
            stop(gettextf("number of offsets (%d) should %d (number of observations)",
                          loff, n), domain = "R-lme4")
        }
    }
    rho$offset <- offset
    if (muetastart) {
        rho$etastart <- unname(as.numeric(model.extract(fr, "etastart")))
        rho$mustart <- unname(as.numeric(model.extract(fr, "mustart")))
    }
}

##' Install derived matrices in the environment.
##'
##' @param rho an lmer environment
##' @fixme Modify for glmer, nlmer, etc. where products are not stored?
derived_mats <- function(rho) {
    stopifnot(is.environment(rho),
              is(X <- rho$X, "Matrix"),
              is(Zt <- rho$Zt, "Matrix"))
    n <- nrow(X)
    p <- ncol(X)
    off <- rho$offset
### The effective y vector has the offset, if present, subtracted
    yy <- rho$y - if (length(off)) off else 0
    stopifnot(ncol(Zt) == n, length(yy) == n)

    ## Zty is stored as a Matrix, not a vector because of a
    ## peculiarity in crossprod(Lambda, Zty) if Lambda is zero
    rho$Zty = Zt %*% yy
    if (p == 0) {
        zz <- c(0L,0L)
        qz <- c(nrow(Zt), 0L)
        rho$Xty <- numeric(0)
        rho$XtX <- new("dpoMatrix", Dim = zz, uplo = "U")
        rho$RX <- new("Cholesky", Dim = zz, uplo = "U", diag = "N")
        rho$ZtX <- new("dgeMatrix", Dim = qz, Dimnames = Zt@Dimnames)
        rho$RZX <- new("dgeMatrix", Dim = qz, Dimnames = Zt@Dimnames)
    } else {
        ## Create crossproduct and check for full column rank
        rho$Xty <- unname(as.vector(crossprod(X, yy)))
        rho$ZtX <- Zt %*% X
        rho$RZX <- solve(rho$L, solve(rho$L, crossprod(rho$Lambda,
                                                       rho$ZtX),
                                      sys = "P"), sys = "L")
        if (is(rho$XtX <- crossprod(X), "dsparseMatrix")) {
            rho$RX <- Cholesky(rho$XtX - crossprod(rho$RZX), LDL = FALSE)
        } else {
            rho$RX <- chol(rho$XtX)
        }
    }
}

##' Install various objects, including Lambda, Lind, Zt and Ut in
##' environment rho based on the list bars of random effects terms and
##' model frame fr.
##'
##' @param bars a list of parsed random-effects terms
##' @param fr a model frame in which to evaluate these terms
##' @param rho an environment that is modified by this function
##' @param checknl logical - should the number of levels be checked
##'
##' @return NULL - the side effect of the function is to modify rho
makeZt <- function(bars, fr, rho, checknl = TRUE) {
    if (!length(bars))
        stop("No random effects terms specified in formula")
    stopifnot(is.list(bars), all(sapply(bars, is.language)),
              inherits(fr, "data.frame"), is.environment(rho))
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
    n <- nrow(fr)
    if(checknl && any(nl >= n))
	stop("Number of levels of a grouping factor for the random effects\n",
	     "must be less than the number of observations")
    ## order terms stably by decreasing number of levels in the factor
    if (any(diff(nl)) > 0) {
        ord <- rev(order(nl))
        blist <- blist[ord]
        nl <- nl[ord]
    }
    rho$Zt <- do.call(rBind, lapply(blist, "[[", "sm"))
    q <- nrow(rho$Zt)

    ## Create and install Lambda, Lind, etc.  This must be done after
    ## any potential reordering of the terms.
    rho$cnms <- lapply(blist, "[[", "cnms")
    nc <- sapply(rho$cnms, length)      # no. of columns per term
    rho$ncTrms <- nc
    nth <- as.integer((nc * (nc+1))/2)  # no. of parameters per term
    nb <- nc * nl                       # no. of random effects per term
    stopifnot(sum(nb) == q)

    rho$u <- numeric(q)
    rho$theta <- numeric(sum(nth))
    boff <- cumsum(c(0L, nb))           # offsets into b
    thoff <- cumsum(c(0L, nth))         # offsets into theta
    rho$Lambda <-
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
    rho$Lind <- as.integer(rho$Lambda@x)
    if (rho$diagonalLambda <- all(nc == 1))
        rho$Lambda <- Diagonal(x = rho$Lambda@x)
    rho$Ut <- crossprod(rho$Lambda, rho$Zt)

    lower <- -Inf * (rho$theta + 1)
    lower[unique(diag(rho$Lambda))] <- 0
    rho$lower <- lower
    rho$theta[] <- is.finite(lower)
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
    rho$flist <- fl

    NULL
} ## {makeZt}

env2lmer <- function(from) {
    ## ---- now "pack" everything into  S4  object --------
    ee <- environment(from@getPars)
    mel <- as.list(ee)                  # <<- via fast internal code
    xnms <- c("verbose", "compDev", "ldRX2", "ldL2", "pwrss",
              "lower", "diagonalLambda",
              "nmp", "p", "weights", "nobs", "REML")
    ## xtrs <- mel[i.xnms <- names(mel) %in% xnms]
    dims <- with(mel[xnms],
                 c(n = nobs, p = p, nmp = nmp, q = nrow(ee$Zt),
                   np = length(ee$theta), # number of parameters in optimization
                   ## logical ones (coerced to 0/1 integer):
                   useSc = TRUE, REML = REML,
                   compDev = compDev, diagonalLambda = diagonalLambda,
                   verb = verbose))
    ## TODO(?): replace devcomp() by direct and the "other" xnms above
    devC <- devcomp(from)
    ## $ cmp = c("ldL2", "ldRX2", "pwrss", "deviance", "REML")
    ## $dims = c("n", "p", "nmp", "q", "useSc")
    mel <- mel[is.na(match(names(mel), xnms))] # == mel[!i.xnms]
    ## enquote(.) the call {TODO: use enquote(.) from R 2.12.0 on}:
    mel[["call"]] <- as.call(list(as.name("quote"), mel[["call"]]))
    ## return  new("mer", .....) :
    do.call(new, c(list("lmer", env = from, dims = dims, devcomp = devC),
                   mel))
}

setAs("lmerenv", "lmer", env2lmer)

lmer <-
    function(formula, data, REML = TRUE, sparseX = FALSE,
             control = list(), start = NULL, verbose = 0, doFit = TRUE,
             compDev = TRUE, optimizer = c("nlminb", "bobyqa", "optimize"),
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
                                        # environment for deviance evaluation
    rho <- new.env(parent = environment(formula))
    rho$REML <- as.logical(REML)[1]
    ## evaluate and install the model frame :
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula)         # substituted "|" by "+" -
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    check_y_weights(fr <- eval(mf, parent.frame()), rho)
    rho$y <- unname(as.double(rho$y))

    fe.form <- formula           # evaluate fixed-effects model matrix
    nb <- nobars(formula[[3]])   # fixed-effects terms only
    if (is.null(nb)) nb <- 1
    fe.form[[3]] <- nb
    X <- if (sparseX) {
        sparse.model.matrix(fe.form, fr, contrasts)
    } else as( model.matrix(fe.form, fr, contrasts), "dgeMatrix")
    rownames(X) <- NULL
    rho$X <- X
    rho$sparseX <- sparseX

    rho$p <- p <- ncol(X)
    stopifnot((rho$nmp <- rho$nobs - p) > 0) # nmp := n m[inus] p
    rho$beta <- numeric(p) ## without names() [speed]

    makeZt(findbars(formula[[3]]), fr, rho)
    rho$L <- Cholesky(tcrossprod(rho$Zt), LDL = FALSE, Imult = 1)

    rho$Lambda@x[] <- rho$theta[rho$Lind]
    derived_mats(rho)

    rho$mu <- numeric(rho$nobs)
    rho$gamma <- numeric(rho$nobs)
    rho$pwrss <- numeric(1)
    rho$ldL2 <- numeric(1)
    rho$ldRX2 <- numeric(1)
    rho$sqrtrWt <- sqrt(rho$weights)
    rho$sqrtXWt <- sqrt(rho$weights)    # to ensure a distinct copy
    ##  sqrtXWt  :=  sqrt of model matrix row weights

    rho$compDev <- compDev
    rho$call <- mc
### The meaning of verbose has changed.
###    verbose == 0 => no iteration output
###    verbose == 1 => iteration output from optimizer but not PIRLS
###    verbose >= 2 => iteration output from optimizer and PIRLS
    rho$verbose <- as.integer(verbose)[1]

    sP <- function(x) {
### FIXME: weights are not yet incorporated (needed here?)
        if (compDev) { ## default: be as fast as possible:
	    ## NOTE: currently need .Call("<name>", ... PACKAGE=*) as this is
	    ##	     called from too low-level (optimizers)
            .Call("lmer_deviance", parent.env(environment()), x,
                  PACKAGE = "lme4a")
        } else {
            stopifnot(length(x) == length(theta))
            ## new theta := x
            theta <<- as.numeric(x)
            Lambda@x[] <<- theta[Lind]  # update Lambda
            Ut <<- crossprod(Lambda, Zt)
            Matrix:::destructive_Chol_update(L, Ut, Imult = 1)
            cu <- solve(L, solve(L, crossprod(Lambda, Zty), sys = "P"),
                        sys = "L")
            RZX <<- solve(L, solve(L, crossprod(Lambda, ZtX), sys = "P"),
                          sys = "L")
            RX <<- if (sparseX) update(RX, XtX - crossprod(RZX)) else
            chol(XtX - crossprod(RZX))
            beta[] <<- if (sparseX) solve(RX, Xty - crossprod(RZX, cu))@x
            else solve(RX, solve(t(RX), Xty - crossprod(RZX, cu)))@x
            u[] <<- solve(L, solve(L, cu - RZX %*% beta, sys = "Lt"),
                          sys = "Pt")@x
            gamma[] <<- (if(length(offset)) offset else 0) +
                (crossprod(Ut, u) + X %*% beta)@x
            mu[] <<- gamma
            pwrss <<- sum(c(y - mu, u)^2) # penalized residual sum of squares
            ldL2[] <<- 2 * determinant(L)$mod
	    if (REML) {
		ldRX2[] <<- 2 * determinant(RX)$mod
		ldL2 + ldRX2 + nmp * (1 + log(2 * pi * pwrss/nmp ))
	    }
	    else
		ldL2	    + nobs * (1 + log(2 * pi * pwrss/nobs))
        }
    }
    gP <- function() theta
    gB <- function() cbind(lower = lower,
                           upper = rep.int(Inf, length(theta)))
    environment(sP) <- environment(gP) <- environment(gB) <- rho
    me <- new("lmerenv", setPars = sP, getPars = gP, getBounds = gB)
    if (doFit) {                        # perform the optimization
        opt <- match.arg(optimizer)
        if (verbose) {
            if (opt == "nlminb") control$trace <- 1
            if (opt == "bobyqa") control$iprint <- 2
            if (opt == "optimize")
                warning("verbose argument ignored with optimize")
        }
        switch(opt,
               bobyqa = bobyqa(me, control = control),
               nlminb = nlminb(me, control = control),
               optimize = if(is.numeric(tl <- control$tol[1]) & tl >= 0)
               optimize(me, tol=tl) else optimize(me))

        env2lmer(me)
    }
    else ## the merenv:
        me
} ## lmer()

glmer <-
function(formula, data, family = gaussian, sparseX = FALSE,
         compDev = TRUE, control = list(), start = NULL,
         verbose = 0, doFit = TRUE, subset, weights, na.action,
         offset, contrasts = NULL, nAGQ = 1, mustart, etastart, ...)
{
    mf <- mc <- match.call()
    if (missing(family)) { ## divert using lmer()
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

    ## environment for deviance evaluation
    rho <- new.env(parent = environment(formula))
    ## evaluate and install the model frame
    m <- match(c("data", "subset", "weights", "na.action", "offset",
                 "mustart", "etastart"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula)
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    check_y_weights(fr <- eval(mf, parent.frame()), rho, muetastart = TRUE)

    ## evaluate fixed-effects model matrix
    fe.form <- formula
    nb <- nobars(formula[[3]])   # fixed-effects terms only
    if (is.null(nb)) nb <- 1
    fe.form[[3]] <- nb
    X <- if (sparseX) {
        sparse.model.matrix(fe.form, fr, contrasts)
    } else as(model.matrix(fe.form, fr, contrasts), "dgeMatrix")
    rownames(X) <- NULL
    rho$X <- X
    rho$sparseX <- sparseX

    p <- ncol(X)
    rho$beta <- numeric(p) ## without names() [speed]

### FIXME: Should decide on the form of the start argument and how it is used
    rho$start <- numeric(p)            # needed for family$initialize

    ## evaluate, check and initialize family
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if(is.function(family)) family <- family()
    rho$family <- family
    ## use scale, whenever not  binomial or poisson :
    rho$useSc <- !match(family$family, c("binomial", "poisson"), nomatch = 0L)

    n <- rho$nobs
    if (!length(rho$weights)) rho$weights <- rep.int(1, n)
    eval(family$initialize, rho)
    rho$weights <- unname(rho$weights)
    mustart <- unname(as.numeric(rho$mustart))
    etastart <- unname(as.numeric(rho$etastart))
    rho$etastart <- rho$mustart <- NULL
    stopifnot(length(etastart) == n || length(mustart) == n)
    if (length(etastart) == n) {
        stopifnot(family$valideta(etastart))
        mustart <- family$linkinv(etastart)
    } else {
        stopifnot(family$validmu(mustart))
        etastart <- family$linkfun(mustart)
    }
    rho$gamma <- etastart
    rho$mu <- mustart
    rho$var <- family$variance(mustart)
    rho$muEta <- family$mu.eta(etastart)
    rho$call <- mc

    ## enforce modes
    rho$y <- unname(as.double(rho$y))   # must be done after initialize
    rho$pwrss <- numeric(1)
    rho$ldL2 <- numeric(1)
    rho$sqrtrWt <- numeric(n)
    rho$sqrtXWt <- numeric(n)
    rho$wtres <- numeric(n)

    checknl <- !(family$family %in% c("binomial", "poisson"))
    makeZt(findbars(formula[[3]]), fr, rho, checknl)
    rho$L <- Cholesky(tcrossprod(rho$Zt), LDL = FALSE, Imult = 1)
    derived_mats(rho)
    rho$compDev <- compDev
    rho$verbose <- as.integer(verbose)[1]
    rho$u[] <- 0                        # for the IRLS step
    rho$lower <- unname(c(rho$lower, rep(-Inf, p)))
    .Call(glmer_IRLS, rho)              # optimize the fixed-effect model

    ## This is the one that will be minimized
    sP <- function(x) {
        rho <- parent.env(environment())
        m <- length(rho$theta)
        sat <- seq_len(m)
        stopifnot(is.numeric(x), length(x) == m + length(beta))
        .Call("merenv_update_Lambda_Ut", rho, x[sat], PACKAGE = "lme4a")
        rho$beta[] <- x[-sat]
        ## return value:
        return(rho$deviance <- .Call("glmer_PIRLS", rho, PACKAGE = "lme4a"))
    }
    gP <- function() c(theta, beta)
    gB <- function() cbind(lower = lower,
                           upper = rep.int(Inf, length(lower)))
    environment(sP) <- environment(gP) <- environment(gB) <- rho
    me <- new("glmerenv", setPars = sP, getPars = gP, getBounds = gB)
    if (doFit) {                        # perform the optimization
        if (verbose) control$trace <- 1
        nlminb(me@getPars(), me@setPars,
               lower = me@getBounds()[,"lower"], control = control)
        .Call(glmer_update_RX, rho)
    }
    me
} ## glmer()

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
nlmer <- function(formula, data, family = gaussian, start = NULL,
                   verbose = 0, nAGQ = 1, doFit = TRUE, subset,
                   weights, na.action, mustart, etastart,
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
                                        # Really do need to check twice
    formula <- as.formula(formula)
    if (length(formula) < 3) stop("formula must be a 3-part formula")
    nlform <- as.formula(formula[[2]])
    if (length(nlform) < 3)
        stop("formula must be a 3-part formula")
    nlmod <- as.call(nlform[[3]])
                                        # check for parameter names in start
    if (is.numeric(start)) start <- list(nlpars = start)
    stopifnot((s <- length(pnames <- names(start$nlpars))) > 0,
              is.numeric(start$nlpars))
    if (!all(pnames %in% (anms <- all.vars(nlmod))))
        stop("not all parameter names are used in the nonlinear model expression")
    fr.form <- nlform
## FIXME: This should be changed to use subbars
    fr.form[[3]] <-         # the frame formula includes all variables
        parse(text = paste(setdiff(all.vars(formula), pnames),
                         collapse = ' + '))[[1]]
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    rho <- new.env(parent = environment(formula))
    fr <- eval(mf, parent.frame())
    check_y_weights(fr, rho)
    attr(rho$frame, "terms") <- NULL
    for (nm in pnames) fr[[nm]] <- start$nlpars[[nm]]
    n <- nrow(fr)

    ## create nlenv and check the evaluation of the nonlinear model function
    rho$nlenv <- new.env()  # want it to inherit from this environment (or formula env)
    lapply(all.vars(nlmod), function(nm) assign(nm, fr[[nm]], envir = rho$nlenv))
    rho$nlmodel <- nlmod
                                        # evaluate and check the family
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if(is.function(family)) family <- family()
    rho$family <- family

    ## these allow starting values to be expressed in terms of other vars.
    rho$mustart <- model.extract(mf, "mustart")
    rho$etastart <- model.extract(mf, "etastart")

    eval(family$initialize, rho)

                                        # enforce modes on some vectors
    rho$y <- unname(as.double(rho$y))
    rho$mustart <- unname(as.double(rho$mustart))
    rho$etastart <- unname(as.double(rho$etastart))
    if (exists("n", envir = rho))
        rho$nobs <- as.double(rho$nobs)

    eta <- eval(rho$nlmodel, rho$nlenv)
    if (is.null(rho$etaGamma <- attr(eta, "gradient")))
        stop("The nonlinear model in nlmer must return a gradient attribute")
    rho$eta <- numeric(length(eta))
    rho$eta[] <- eta

    ## build the extended frame for evaluation of X and Zt
    fr <- do.call(rbind, lapply(1:s, function(i) fr)) # rbind s copies of the frame
    for (nm in pnames) # convert these variables in fr to indicators
        fr[[nm]] <- as.numeric(rep(nm == pnames, each = n))
    fe.form <- nlform # modify formula to suppress intercept (Is this a good idea?)
    fe.form[[3]] <- substitute(0 + bar, list(bar = nobars(formula[[3]])))
    rho$X <- model.matrix(fe.form, fr, contrasts)
    rownames(rho$X) <- NULL
    p <- ncol(rho$X)
    if ((qrX <- qr(rho$X))$rank < p)
        stop(gettextf("rank of X = %d < ncol(X) = %d", qrX$rank, p))
    rho$start <- numeric(p)  # must be careful that these are distinct
    rho$start[] <- rho$beta <- qr.coef(qrX, unlist(lapply(pnames, get, envir = rho$nlenv)))
    rho$RX <- qr.R(qrX)
#    eb <- evalbars(formula, fr, contrasts, TRUE) # flist, trms, nest
#    rho$dims["nest"] <- eb$nest
#    rho$flist <- eb$flist
#    lmerFactorList(eb$trms, rho)

    q <- length(rho$u)
    rho$u0 <- numeric(q)
                                        # evaluate the control argument
#    control <- do.call(lmerControl, as.list(control))
    if (!missing(verbose)) control$trace <- as.integer(verbose[1])
    rho$dims["verb"] <- control$trace
    control$trace <- abs(control$trace) # negative values give PIRLS output
    rho$control <- control
    rho$mc <- mc                        # store the matched call

    rho$bds <- getBounds(rho)
    theta0 <- getPars(rho)
    if (!is.null(start$theta)) {
        stopifnot(length(st <- as.double(start$theta)) == length(theta0))
        setPars(rho, st)
    } else {
        setPars(rho, theta0) # one evaluation to check structure
    }
    rho$beta0[] <- rho$beta
    rho$u0[] <- rho$u
    if (!doFit) return(rho)

#    merFinalize(rho)
}

setMethod("formula", "merenv", function(x, ...) formula(env(x)$call, ...))
setMethod("fixef", "merenv", function(object, ...) {
    e <- env(object)
    structure(e$beta, names = dimnames(e$X)[[2]]) })

setMethod("formula", "mer", function(x, ...) formula(x@call, ...))
setMethod("fixef", "mer", function(object, ...)
	  structure(object@beta, names = dimnames(object@X)[[2]]))

setMethod("formula", "lmerTrms", function(x, ...) formula(x@call, ...))
setMethod("fixef", "deFeMod", function(object, ...)
	  structure(object@beta, names = dimnames(object@X)[[2]]))
setMethod("fixef", "lmerTrms", function(object, ...) fixef(object@fe))


if(FALSE) {## These are not used, rather the methods below
setMethod("ranef", "merenv", function(object, ...)
	  with(env(object), (Lambda %*% u)@x))
setMethod("ranef", "mer", function(object, ...)
	  (object@Lambda %*% object@u)@x)
}

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

setMethod("ranef", signature(object = "merenvtrms"),
          function(object, postVar = FALSE, drop = FALSE,
                   whichel = names(ans), ...)
      {
          rho <- env(object)
          ## evaluate the list of matrices
          levs <- lapply(fl <- rho$flist, levels)
          asgn <- attr(fl, "assign")
          cnms <- rho$cnms
          nc <- rho$ncTrms
          nb <- nc * (nl <- unlist(lapply(levs, length))[asgn])
          nbseq <- rep.int(seq_along(nb), nb)
          ml <- split(as.vector(rho$Lambda %*% rho$u), nbseq)
          for (i in seq_along(ml))
              ml[[i]] <- matrix(ml[[i]], nc = nc[i], byrow = TRUE,
                                dimnames = list(NULL, cnms[[i]]))
          ## produce a list of data frames corresponding to
          ## factors, not terms
          ans <- lapply(seq_along(fl),
                        function(i)
                        data.frame(do.call(cbind, ml[asgn == i]),
                                   row.names = levs[[i]],
                                   check.names = FALSE))
          names(ans) <- names(fl)

          ## Process whichel
          stopifnot(is(whichel, "character"))
          whchL <- names(ans) %in% whichel
          ans <- ans[whchL]

          if (postVar) {
              vv <- .Call(merenvtrms_condVar, rho, sigma(object))
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
          ans
      })

## This is a BIG CHEAT [ as the above uses C++ ....]  __ FIXME __
## is *wrong* when the 'env' has been changed
## and not really efficient {have already copied everything into slots}:
setMethod("ranef", signature(object = "mer"),
	  function(object, postVar = FALSE, drop = FALSE, ...)
	  ranef(object@env, postVar=postVar, drop=drop, ...))
##              ----------

dotplot.ranef.mer <- function(x, data, ...)
{
    prepanel.ci <- function(x, y, se, subscripts, ...) {
        if (is.null(se)) return(list())
        x <- as.numeric(x)
        hw <- 1.96 * as.numeric(se[subscripts])
        list(xlim = range(x - hw, x + hw, finite = TRUE))
    }
    panel.ci <- function(x, y, se, subscripts, pch = 16,
                         horizontal = TRUE, col = dot.symbol$col,
                         lty = dot.line$lty, lwd = dot.line$lwd,
                         col.line = dot.line$col, levels.fos = unique(y),
                         groups = NULL, ...)
    {
        x <- as.numeric(x)
        y <- as.numeric(y)
        dot.line <- trellis.par.get("dot.line")
        dot.symbol <- trellis.par.get("dot.symbol")
        sup.symbol <- trellis.par.get("superpose.symbol")
        panel.abline(h = levels.fos, col = col.line, lty = lty, lwd = lwd)
        panel.abline(v = 0, col = col.line, lty = lty, lwd = lwd)
        if (!is.null(se)) {
            se <- as.numeric(se[subscripts])
            panel.segments( x - 1.96 * se, y, x + 1.96 * se, y, col = 'black')
        }
        panel.xyplot(x, y, pch = pch, ...)
    }
    f <- function(x, ...) {
        ss <- stack(x)
        ss$ind <- factor(as.character(ss$ind), levels = colnames(x))
        ss$.nn <- rep.int(reorder(factor(rownames(x)), x[[1]]), ncol(x))
        se <- NULL
        if (!is.null(pv <- attr(x, "postVar")))
            se <- unlist(lapply(1:(dim(pv)[1]), function(i) sqrt(pv[i, i, ])))
        dotplot(.nn ~ values | ind, ss, se = se,
                prepanel = prepanel.ci, panel = panel.ci,
                xlab = NULL, ...)
    }
    lapply(x, f, ...)
}

plot.ranef.mer <- function(x, y, ...)
{
    lapply(x, function(x) {
        cn <- lapply(colnames(x), as.name)
        switch(min(ncol(x), 3),
               qqmath(eval(substitute(~ x, list(x = cn[[1]]))), x, ...),
               xyplot(eval(substitute(y ~ x,
                                      list(y = cn[[1]],
                                           x = cn[[2]]))), x, ...),
               splom(~ x, ...))
    })
}

qqmath.ranef.mer <- function(x, data, ...)
{
    prepanel.ci <- function(x, y, se, subscripts, ...) {
        y <- as.numeric(y)
        se <- as.numeric(se[subscripts])
        hw <- 1.96 * se
        list(ylim = range(y - hw, y + hw, finite = TRUE))
    }
    panel.ci <- function(x, y, se, subscripts, pch = 16, ...)  {
        panel.grid(h = -1,v = -1)
        panel.abline(h = 0)
        x <- as.numeric(x)
        y <- as.numeric(y)
        se <- as.numeric(se[subscripts])
        ly <- y - 1.96 * se
        uy <- y + 1.96 * se
        panel.segments(x, y - 1.96*se, x, y + 1.96 * se,
                       col = 'black')
        panel.xyplot(x, y, pch = pch, ...)
    }
    f <- function(x) {
        if (!is.null(pv <- attr(x, "postVar"))) {
            cols <- 1:(dim(pv)[1])
            se <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
            nr <- nrow(x)
            nc <- ncol(x)
            ord <- unlist(lapply(x, order)) +
                rep((0:(nc - 1)) * nr, each = nr)
            rr <- 1:nr
            ind <- gl(ncol(x), nrow(x), labels = names(x))
            xyplot(unlist(x)[ord] ~
                   rep(qnorm((rr - 0.5)/nr), ncol(x)) | ind[ord],
                   se = se[ord], prepanel = prepanel.ci, panel = panel.ci,
                   scales = list(y = list(relation = "free")),
                   xlab = "Standard normal quantiles",
                   ylab = NULL, aspect = 1, ...)
        } else {
            qqmath(~values|ind, stack(x),
                   scales = list(y = list(relation = "free")),
                   xlab = "Standard normal quantiles",
                   ylab = NULL, ...)
        }
    }
    lapply(x, f)
}

plot.coef.mer <- function(x, y, ...)
{
    ## remove non-varying columns from frames
    reduced <- lapply(x, function(el)
                      el[, !sapply(el, function(cc) all(cc == cc[1]))])
    plot.ranef.mer(reduced, ...)
}

dotplot.coef.mer <- function(x, data, ...) {
    mc <- match.call()
    mc[[1]] <- as.name("dotplot.ranef.mer")
    eval(mc)
}



.devc.lmer <- function(nobs, nmp, ldL2, ldRX2, pwrss) {
    c(deviance = ldL2        + nobs * (1 + log(2 * pi * pwrss/nobs)),
      REML     = ldL2 + ldRX2 + nmp * (1 + log(2 * pi * pwrss/nmp)))
}

setMethod("devcomp", "lmerenv", function(x, ...)
	  with(env(x),# somewhat "ugly" (need 'lme4a:::' below ...)
	       list(cmp =
		    c(ldL2 = ldL2, ldRX2 = ldRX2, pwrss = pwrss,
                      lme4a:::.devc.lmer(nobs, nmp, ldL2, ldRX2, pwrss)),
		    dims =
		    c(n = nobs, p = length(beta), nmp = nmp, q = nrow(Zt),
		      useSc = TRUE))))

## Cheap for now; not sure if the slot should be kept at all:
setMethod("devcomp", "mer", function(x, ...) x@devcomp)

setMethod("devcomp", "lmerTrms", function(x, ...)
      {
	  n <- nrow  (x@fe @ X)
	  p <- length(x@fe @ beta)
	  ldRX2 <-    x@fe @ ldRX2
	  ldL2 <- x@re @ ldL2
	  u    <- x@re @ u
	  pwrss <- x@resp@wrss + sum(u*u)
	  list(cmp =
	       c(ldL2 = ldL2, ldRX2 = ldRX2, pwrss = pwrss,
		 .devc.lmer(n, n-p, ldL2, ldRX2, pwrss)),
	       dims = c(n = n, p = p, nmp = n-p, q = length(u), useSc = TRUE))
      })


setMethod("devcomp", "glmerenv",
	  function(x, ...)
	  with(env(x),
	       list(cmp = c(ldL2 = ldL2, deviance = deviance),
		    dims = c(n = nobs, p = length(beta), q = nrow(Zt),
			     useSc = useSc))))


if(FALSE) ## I don't think this is used
setMethod("env", "mer", function(x) env(x@env))

setMethod("getL", "merenv", function(x) env(x)$L)
setMethod("getL", "mer", function(x) x@L)
setMethod("getL", "reModule", function(x) x@L)
setMethod("getL", "lmer2", function(x) x@re@L)


setMethod("sigma", signature(object = "merenv"),
	  function (object, ...) {
	      dc <- devcomp(object)
	      if(dc$dims[["useSc"]])
		  sqrt(dc$cmp[["pwrss"]]/
		       dc$dims[[if (env(object)$REML) "nmp" else "n"]])
	      else 1
	  })
setMethod("sigma", signature(object = "mer"), function (object, ...)
      {
	  dm <- object@dims
	  if(dm[["useSc"]])
	      sqrt(object@devcomp$cmp[["pwrss"]]/
		   dm[[if(dm[["REML"]]) "nmp" else "n"]])
	  else 1
      })

setMethod("sigma", signature(object = "lmerTrms"), function (object, ...)
      {
	  dc <- devcomp(object)
	  dm <- dc$dims
## in general {maybe use, making this into aux.function?
## 	  if(dm[["useSc"]])
## 	      sqrt(dc$cmp[["pwrss"]]/
## 		   dm[[if(dc$cmp[["REML"]]) "nmp" else "n"]])
## 	  else 1
          sqrt(dc$cmp[["pwrss"]]/ dm[[if(object@REML) "nmp" else "n"]])
      })


## Obtain the ML fit of the model parameters
MLfit <- function(x) {
    stopifnot(is(x, "lmerenv"))
    if (!env(x)$REML) return(x)
    update(x, REML = FALSE)
}

setMethod("isREML", "lmerenv",	function(x) isTRUE(env(x)$REML))
setMethod("isREML", "lmer",	function(x) as.logical(x@dims[["REML"]]))
setMethod("isREML", "merenv",	function(x) FALSE)
setMethod("isREML", "mer",	function(x) FALSE)

setMethod("getCall", "merenv",	function(x) env(x)$call)
setMethod("getCall", "mer",	function(x) x@call)

##' <description>
##'
##' <details>
##' @title anova() for both  "lmerenv" and "lmer" fitted models
##' @param object an "lmerenv" or "lmer" - fitted model
##' @param ...  further such objects
##' @return an "anova" data frame; the traditional (S3) result of anova()
anovaLmer <- function(object, ...) {
    mCall <- match.call(expand.dots = TRUE)
    dots <- list(...)
    .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))
    modp <- {
	as.logical(.sapply(dots, is, "lmerenv")) |
	as.logical(.sapply(dots, is, "lmer")) |
	as.logical(.sapply(dots, is, "lm")) }
    if (any(modp)) {			# multiple models - form table
	opts <- dots[!modp]
	mods <- c(list(object), dots[modp])
	## model names
	mNms <- .sapply(as.list(mCall)[c(FALSE, TRUE, modp)], deparse)
	names(mods) <- sub("@env$", '', mNms) # <- hack
	mods <- lapply(mods, function(x) ## get ML estimates if necessary
		   {
		       if (isREML(x)) return(x)
		       update(x, REML = FALSE)
		   })

	llks <- lapply(mods, logLik)
        ii <- order(Df <- .sapply(llks, attr, "df"))
	mods <- mods[ii]
	llks <- llks[ii]
	Df   <- Df  [ii]
	calls <- lapply(mods, getCall)
	data <- lapply(calls, "[[", "data")
	if (any(data != data[[1]]))
	    stop("all models must be fit to the same data object")
	header <- paste("Data:", data[[1]])
	subset <- lapply(calls, "[[", "subset")
	if (any(subset != subset[[1]]))
	    stop("all models must use the same subset")
	if (!is.null(subset[[1]]))
	    header <-
		c(header, paste("Subset", deparse(subset[[1]]),
				sep = ": "))
	llk <- unlist(llks)
	chisq <- 2 * pmax(0, c(NA, diff(llk)))
	dfChisq <- c(NA, diff(Df))
	val <- data.frame(Df = Df,
			  AIC = .sapply(llks, AIC),
			  BIC = .sapply(llks, BIC),
			  logLik = llk,
			  "Chisq" = chisq,
			  "Chi Df" = dfChisq,
			  "Pr(>Chisq)" = pchisq(chisq, dfChisq, lower = FALSE),
			  row.names = names(mods), check.names = FALSE)
	class(val) <- c("anova", class(val))
	attr(val, "heading") <-
	    c(header, "Models:",
	      paste(rep(names(mods), times = unlist(lapply(lapply(lapply(calls,
				     "[[", "formula"), deparse), length))),
		    unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
		    sep = ": "))
	return(val)
    }
    else { ## ------ single model ---------------------
	dc <- devcomp(object)
	p <- dc$dims[["p"]]
	ss <- fixef(object)
	## FIXME :
	stop("assign attribute not currently available")
	asgn <- attr(object@X, "assign")
	terms <- terms(object)
	nmeffects <- attr(terms, "term.labels")
	if ("(Intercept)" %in% names(ss))
	    nmeffects <- c("(Intercept)", nmeffects)
	ss <- unlist(lapply(split(ss, asgn), sum))
	df <- unlist(lapply(split(asgn,	 asgn), length))
	## dfr <- unlist(lapply(split(dfr, asgn), function(x) x[1]))
	ms <- ss/df
	f <- ms/(sigma(object)^2)
	## P <- pf(f, df, dfr, lower.tail = FALSE)
	## table <- data.frame(df, ss, ms, dfr, f, P)
	table <- data.frame(df, ss, ms, f)
	dimnames(table) <-
	    list(nmeffects,
		 ## c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
		 c("Df", "Sum Sq", "Mean Sq", "F value"))
	if ("(Intercept)" %in% nmeffects)
	    table <- table[-match("(Intercept)", nmeffects), ]
	attr(table, "heading") <- "Analysis of Variance Table"
	class(table) <- c("anova", "data.frame")
	table
    }
}

setMethod("anova", signature(object = "lmerenv"), anovaLmer)

setMethod("anova", signature(object = "lmer"), anovaLmer)

##' <description>
##'
##' <details>
##' @title vcov(): Extract conditional covariance matrix of fixed effects
mkVcov <- function(object, RX, nmsX, correlation = TRUE, ...) {
    V <- sigma(object)^2 * chol2inv(RX)
    if(is.null(rr <- tryCatch(as(V, "dpoMatrix"),
                              error = function(e) NULL)))
        stop("Computed variance-covariance matrix is not positive definite")
    dimnames(rr) <- list(nmsX, nmsX)
    if(correlation)
        rr@factors$correlation <- as(rr, "corMatrix")
    rr
}
##' Extract the conditional variance-covariance matrix of the fixed effects
##' @param object
##' @param correlation
##' @param ...
setMethod("vcov", signature(object = "merenv"),
	  function(object, correlation = TRUE, ...)
      {
	  en <- env(object)
	  mkVcov(object, RX = en$RX, nmsX = colnames(en$X),
		 correlation=correlation, ...)
      })

setMethod("vcov", signature(object = "mer"),
	  function(object, correlation = TRUE, ...)
	  mkVcov(object, RX = object@RX, nmsX = colnames(object@X),
		 correlation=correlation, ...))


mkVarCorr <- function(sc, cnms, nc, theta, flist) {
    ncseq <- seq_along(nc)
    thl <- split(theta, rep.int(ncseq, (nc * (nc + 1))/2))
    ans <- lapply(ncseq, function(i)
              {
                  mm <- diag(nrow = nc[i])
                  mm[lower.tri(mm, diag = TRUE)] <- thl[[i]]
                  rownames(mm) <- cnms[[i]]
                  val <- tcrossprod(sc * mm) # variance-covariance
                  stddev <- sqrt(diag(val))
                  correl <- t(val / stddev)/stddev
                  diag(correl) <- 1
                  attr(val, "stddev") <- stddev
                  attr(val, "correlation") <- correl
                  val
              })
    names(ans) <- names(flist)[attr(flist, "assign")]
    attr(ans, "sc") <- sc
    ans
}
##' Create the VarCorr object of variances and covariances
setMethod("VarCorr", signature(x = "merenv"), function(x) {
    rho <- env(x)
    mkVarCorr(sigma(x), cnms=rho$cnms, nc = rho$ncTrms,
	      theta=rho$theta, flist=rho$flist)
})

setMethod("VarCorr", signature(x = "mer"), function(x)
    mkVarCorr(sigma(x), cnms=x@cnms, nc = x@ncTrms,
	      theta=x@theta, flist=x@flist))

setMethod("VarCorr", signature(x = "lmerTrms"), function(x) {
    cnms <- x@trms@cnms
    nc <- sapply(cnms, length)      # no. of columns per term
    mkVarCorr(sigma(x), cnms=cnms, nc = nc,
	      theta=x@re@theta, flist=x@trms@flist)
})


## This is modeled a bit after  print.summary.lm :
## Prints *both*  'mer' and 'merenv' - as it uses summary(x) mainly
printMerenv <- function(x, digits = max(3, getOption("digits") - 3),
                        correlation = NULL, symbolic.cor = FALSE,
                        signif.stars = getOption("show.signif.stars"), ...)
{
    so <- summary(x)
    cat(sprintf("%s ['%s']\n",so$methTitle, class(x)))
    if (!is.null(f <- so$family))
	cat(" Family:", f$family,"\n")
    if (!is.null(cc <- so$call$formula))
	cat("Formula:", deparse(cc),"\n")
    if (!is.null(cc <- so$call$data))
	cat("   Data:", deparse(cc), "\n")
    if (!is.null(cc <- so$call$subset))
	cat(" Subset:", deparse(asOneSidedFormula(cc)[[2]]),"\n")
    ## MM would like:    cat("\n")
    print(so$AICtab, digits = digits)
    cat("\nRandom effects:\n")
    print(formatVC(so$varcor, digits = digits, useScale = so$useScale),
	  quote = FALSE, digits = digits, ...)

    ngrps <- so$ngrps
    cat(sprintf("Number of obs: %d, groups: ", so$devcomp$dims[["n"]]))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
    p <- nrow(so$coefficients)
    if (p > 0) {
	cat("\nFixed effects:\n")
	printCoefmat(so$coefficients, zap.ind = 3, #, tst.ind = 4
		     digits = digits, signif.stars = signif.stars)
	if(!is.logical(correlation)) { # default
	    correlation <- p <= 20
	    if(!correlation) {
		nam <- deparse(substitute(x)) # << TODO: improve if this is called from show()
		cat(sprintf(paste("\nCorrelation matrix not shown, as p = %d > 20.",
				  "Use print(%s, correlation=TRUE)  or",
				  "    vcov(%s)	 if you need it\n", sep="\n"),
			    p, nam, nam))
	    }
	}
	if(correlation) {
	    if(is.null(VC <- so$vcov)) VC <- vcov(x)
	    corF <- VC@factors$correlation
	    if (is.null(corF)) {
		cat("\nCorrelation of Fixed Effets is not available\n")
	    }
	    else {
		p <- ncol(corF)
		if (p > 1) {
		    rn <- rownames(so$coefficients)
		    rns <- abbreviate(rn, minlen=11)
		    cat("\nCorrelation of Fixed Effects:\n")
		    if (is.logical(symbolic.cor) && symbolic.cor) {
			corf <- as(corF, "matrix")
			dimnames(corf) <- list(rns,
					       abbreviate(rn, minlength=1, strict=TRUE))
			print(symnum(corf))
		    }
		    else {
			corf <- matrix(format(round(corF@x, 3), nsmall = 3),
				       nc = p,
                                       dimnames = list(rns, abbreviate(rn, minlen=6)))
			corf[!lower.tri(corf)] <- ""
			print(corf[-1, -p, drop=FALSE], quote = FALSE)
		    }
		}
	    }
	}
    }
    invisible(x)
}## printMerenv()

##' <description>
##' Compute standard errors of fixed effects from an lmer()
##'
##' <details>
##' @title
##' @param object "lmerenv" object,
##' @param RX the Cholesky factor (CHMfactor) L of ...
##' @return numeric vector of length length(fixef(.))
##' @author Doug Bates & Martin Maechler
##' currently *not* exported on purpose
unscaledVar <- function(object, RX = env(object)$RX)
{
    if (is(RX, "Cholesky")) return(diag(chol2inv(RX)))
    stopifnot(is(RX, "CHMfactor"))
    p <- ncol(RX)
    if (p < 1) return(numeric(0))
    p1 <- p - 1L
    ei <- as(c(1, rep.int(0, p1)), "sparseMatrix")
    DI <- function(i) {
	ei@i <- i
	sum(solve(RX, ei, sys = "L")@x^2)
    }
    as.vector(solve(RX, unlist(lapply(0:p1, DI)),
		    system = "Pt"))
}

formatVC <- function(varc, digits = max(3, getOption("digits") - 2),
                     useScale)
### "format()" the 'VarCorr' matrix of the random effects -- for show()ing
{
    sc <- unname(attr(varc, "sc"))
    recorr <- lapply(varc, attr, "correlation")
    reStdDev <- c(lapply(varc, attr, "stddev"), if(useScale) list(Residual = sc))
    reLens <- unlist(c(lapply(reStdDev, length)))
    nr <- sum(reLens)
    reMat <- array('', c(nr, 4),
		   list(rep.int('', nr),
			c("Groups", "Name", "Variance", "Std.Dev.")))
    reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
    reMat[,2] <- c(unlist(lapply(varc, colnames)), if(useScale) "")
    reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
    reMat[,4] <- format(unlist(reStdDev), digits = digits)
    if (any(reLens > 1)) {
	maxlen <- max(reLens)
	corr <-
	    do.call("rBind",
		    lapply(recorr,
			   function(x, maxlen) {
			       x <- as(x, "matrix")
			       cc <- format(round(x, 3), nsmall = 3)
			       cc[!lower.tri(cc)] <- ""
			       nr <- dim(cc)[1]
			       if (nr >= maxlen) return(cc)
			       cbind(cc, matrix("", nr, maxlen-nr))
			   }, maxlen))
	colnames(corr) <- c("Corr", rep.int("", maxlen - 1))
	cbind(reMat, rBind(corr, rep.int("", ncol(corr))))
    } else reMat
}

.summMer <- function(x, rho, devC, flags, varcov = FALSE)
{
    ## Compute  'vcov' only when  'varcov = TRUE', no longer by default

    ## cld <- getClass(cl <- class(object))
    ## isLmer <- extends(cld, "lmerenv")
    ## isGLmer <- extends(cld, "glmerenv")

    ## devC <- devcomp(x)
    ## rho <- env(x)
    isLmer <- flags[["isLmer"]]
    REML <- if(isLmer) flags[["REML"]] else FALSE

    useSc <- as.logical(devC$dims["useSc"])
    fcoef <- fixef(x)
    sig <- if(useSc) sigma(x) else 1.0
    coefs <- cbind("Estimate" = fcoef,
                   "Std. Error" = sig * sqrt(unscaledVar(RX = rho$RX)))
    if (nrow(coefs) > 0) {
        coefs <- cbind(coefs, coefs[,1]/coefs[,2], deparse.level=0)
        colnames(coefs)[3] <- paste(if(useSc) "t" else "z", "value")
    }
    llik <- logLik(x, REML)
    mName <- paste(if(isLmer) "Linear" else if(flags[["isGLmer"]]) "Generalized linear"
                   else "Nonlinear (?)",# TODO: nonlinear / generalized nonlin
                   "mixed model fit by",
                   if (REML) "REML" else "maximum likelihood",
                   if(flags[["isGLmer"]]) "(Laplace Approximation)")
    AICstats <- {
        if (REML) devC$cmp["REML"] # do *not* show likelihood stats here
        else c(AIC = AIC(llik), BIC = BIC(llik), logLik = c(llik),
               deviance = devC$cmp[["deviance"]])
    }
    varcor <- VarCorr(x)

    ## use S3 class for now
    structure(list(methTitle = mName,
                   devcomp = devC, isLmer = isLmer, useScale = useSc,
                   logLik = llik, family = rho$family,
                   ngrps = sapply(rho$flist, function(x) length(levels(x))),
                   coefficients = coefs,
                   sigma = sig,
                   vcov = if(varcov) vcov(x),
                   varcor = varcor, # and use formatVC(.) for printing.
                   AICtab= AICstats,
                   call = rho$call
                   ), class = "summary.merenv")
} ## .summMer()

summaryMerenv <- function(object, varcov = FALSE, ...)
{
    ## Compute  'vcov' only when  'varcov = TRUE', no longer by default

    cld <- getClass(cl <- class(object))
    rho <- env(object)
    .summMer(object, rho=rho, devC = devcomp(object),
	     flags = c(REML=rho$REML,
			isLmer = extends(cld, "lmerenv"),
			isGLmer= extends(cld, "glmerenv")),
	     varcov=varcov)
} ## summaryMerenv()

summaryMer <- function(object, varcov = FALSE, ...)
{
    ## Compute	'vcov' only when  'varcov = TRUE', no longer by default

    cld <- getClass(cl <- class(object))
    isG <- extends(cld, "glmer")
    rho <- sapply(c("RX", if(isG) "family", "flist", "call"),
		  slot, object=object,
		  simplify = FALSE)
    .summMer(object, rho=rho, devC = object@devcomp,
	     flags = c(REML = object@dims[["REML"]],
			isLmer = extends(cld, "lmer"), isGLmer= isG),
	     varcov=varcov)
} ## summaryMer()

summaryMer2 <- function(object, varcov = FALSE, ...)
{
    ## Compute	'vcov' only when  'varcov = TRUE', no longer by default

    cld <- getClass(cl <- class(object))
    ## begin{workaround} -- FIXME (later)
    stopifnot(extends(cld, "lmerTrmsDe"))
    isG <- FALSE ## extends(cld, "glmer")
    isLmer <- TRUE
    ## end{workaround}

    rho <- list(RX = object@fe@RX,
                ## FIXME : "family" !
                flist = object@trms@flist,
                call = object@call)
    devC <- devcomp(object)
    .summMer(object, rho=rho, devC = devC,
	     flags = c(REML = object@REML, isLmer = isLmer, isGLmer= isG),
	     varcov=varcov)
} ## summaryMer2()


setMethod("summary", "lmerenv", summaryMerenv)
## and the *same* here [just for now?]:
setMethod("summary", "glmerenv", summaryMerenv)

setMethod("summary", "mer", summaryMer)

setMethod("summary", "lmerTrms", summaryMer2)


setMethod("model.frame", signature(formula = "merenv"),
	  function(formula, ...) env(formula)$frame)

setMethod("model.matrix", signature(object = "merenv"),
	  function(object, ...) env(object)$X)

setMethod("terms", signature(x = "merenv"),
	  function(x, ...) attr(env(x)$frame, "terms"))

setMethod("model.frame", signature(formula = "mer"),
	  function(formula, ...) formula$frame)
setMethod("model.matrix", signature(object = "mer"),
	  function(object, ...) object$X)
setMethod("terms", signature(x = "mer"),
	  function(x, ...) attr(x$frame, "terms"))

setMethod("deviance", signature(object="lmerenv"),
	  function(object, REML = NULL, ...)
      {
          if (missing(REML) || is.null(REML) || is.na(REML[1]))
              REML <- env(object)$REML
	  devcomp(object)$cmp[[if(REML) "REML" else "deviance"]]
      })
setMethod("deviance", signature(object="lmerTrms"),
	  function(object, REML = NULL, ...)
      {
          if (missing(REML) || is.null(REML) || is.na(REML[1]))
              REML <- object@REML
	  devcomp(object)$cmp[[if(REML) "REML" else "deviance"]]
      })
setMethod("deviance", signature(object="lmer"),
	  function(object, REML = NULL, ...)
      {
	  if (missing(REML) || is.null(REML) || is.na(REML[1]))
	      REML <- object@dims[["REML"]]
	  object@devcomp$cmp[[if(REML) "REML" else "deviance"]]
      })

## The non-lmerenv case:
setMethod("deviance", signature(object="merenv"),
	  function(object, ...) devcomp(object)$cmp[["deviance"]])

## The non-lmer case:
setMethod("deviance", signature(object="mer"),
	  function(object, ...) object@devcomp$cmp[["deviance"]])


mkLogLik <- function(dev, n, p, np, REML = NULL, hasScale) {
    val <- - dev/2
    attr(val, "nall") <- attr(val, "nobs") <- n
    attr(val, "df") <- p + np + hasScale
    if(!is.null(REML))
	attr(val, "REML") <-  as.logical(REML)
    class(val) <- "logLik"
    val
}
##' Extract the log-likelihood or restricted log-likelihood
setMethod("logLik", signature(object="lmerenv"),
	  function(object, REML = NULL, ...)
      {
	  rho <- env(object)
	  if (is.null(REML) || is.na(REML[1]))
	      REML <- rho$REML
	  mkLogLik(deviance(object, REML = REML),
		   n = length(rho$y),
		   p = length(rho$beta), np = length(rho$theta),
		   REML = REML, hasScale = 1L)
      })
setMethod("logLik", signature(object="lmer"),
	  function(object, REML = NULL, ...)
      {
	  dm <- object@dims
	  if (is.null(REML) || is.na(REML[1]))
	      REML <- dm[["REML"]]
	  mkLogLik(deviance(object, REML = REML),
		   n = dm[["n"]], p = dm[["p"]], np = dm[["np"]],
		   REML = REML, hasScale = 1L)
      })
setMethod("sigma", signature(object = "lmerTrms"), function (object, ...)
      {
	  dc <- devcomp(object)
	  dm <- dc$dims
## in general {maybe use, making this into aux.function?
## 	  if(dm[["useSc"]])
## 	      sqrt(dc$cmp[["pwrss"]]/
## 		   dm[[if(dc$cmp[["REML"]]) "nmp" else "n"]])
## 	  else 1
          sqrt(dc$cmp[["pwrss"]]/ dm[[if(object@REML) "nmp" else "n"]])
      })

setMethod("logLik", signature(object="lmerTrms"),
	  function(object, REML = NULL, ...)
      {
	  if (is.null(REML) || is.na(REML[1]))
	      REML <- object@REML
	  mkLogLik(deviance(object, REML = REML),
		   n  = nrow  (object@fe @ X),
		   p  = length(object@fe @ beta),
		   np = length(object@re @ theta),
		   REML = REML, hasScale = 1L)
      })

## The non-lmerenv case:
setMethod("logLik", signature(object="merenv"), function(object, ...)
      {
	  rho <- env(object)
	  stopifnot(is.logical(useSc <- rho$useSc))# required for non-lmerenv
	  mkLogLik(deviance(object),
		   n = length(rho$y),
		   p = length(rho$beta), np = length(rho$theta),
		   hasScale = useSc)
      })
setMethod("logLik", signature(object="mer"),
	  function(object, REML = NULL, ...)
      {
	  dm <- object@dims
	  stopifnot(is.logical(useSc <- dm[["useSc"]]))# required for non-lmer
	  mkLogLik(deviance(object),
		   n = dm[["n"]], p = dm[["p"]], np = dm[["np"]],
		   hasScale = useSc)
      })

setMethod("logLik", signature(object="mer"),
	  function(object, REML = NULL, ...)
      {
	  dm <- object@dims
	  stopifnot(is.logical(useSc <- dm[["useSc"]]))# required for non-lmer
	  mkLogLik(deviance(object),
		   n = dm[["n"]], p = dm[["p"]], np = dm[["np"]],
		   hasScale = useSc)
      })

##' update()  both for "merenv" and "mer" :
updateMer <- function(object, formula., ..., evaluate = TRUE)
{
    if (is.null(call <- getCall(object)))
	stop("object should contain a 'call' component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.))
	call$formula <- update.formula(formula(object), formula.)
    if (length(extras) > 0) {
	existing <- !is.na(match(names(extras), names(call)))
	for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
	if (any(!existing)) {
	    call <- c(as.list(call), extras[!existing])
	    call <- as.call(call)
	}
    }
    if (evaluate)
	eval(call, parent.frame())
    else call
}

setMethod("update", signature(object = "merenv"), updateMer)
setMethod("update", signature(object = "mer"), updateMer)


setMethod("print", "merenv", printMerenv)
setMethod("show", "merenv", function(object) printMerenv(object))

setMethod("print", "mer", printMerenv)
setMethod("show", "mer", function(object) printMerenv(object))

coefMer <- function(object, ...)
{
    if (length(list(...)))
        warning(paste('arguments named "',
                      paste(names(list(...)), collapse = ", "),
                      '" ignored', sep = ''))
    fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
    ref <- ranef(object)
    val <- lapply(ref, function(x)
                  fef[rep.int(1L, nrow(x)),,drop = FALSE])
    for (i in seq(a = val)) {
        refi <- ref[[i]]
        row.names(val[[i]]) <- row.names(refi)
        nmsi <- colnames(refi)
        if (!all(nmsi %in% names(fef)))
            stop("unable to align random and fixed effects")
        for (nm in nmsi) val[[i]][[nm]] <- val[[i]][[nm]] + refi[,nm]
    }
    class(val) <- "coef.mer"
    val
} ##  {coefMer}

setMethod("coef", signature(object = "merenvtrms"), coefMer)
## for now
setMethod("coef", signature(object = "mer"), coefMer)


## For Matrix API change (Oct.2009):
assign("det_CHMfactor.warn", FALSE, envir = Matrix:::.MatrixEnv)
