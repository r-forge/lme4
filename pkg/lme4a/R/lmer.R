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

##' Install derived matrices in the environment.  The y vector is
##' passed separately so that it can incorporate an offset.
##'
##' @param rho an lmer environment
derived_mats <- function(rho) {
    stopifnot(is.environment(rho),
              is(X <- rho$X, "Matrix"),
              is(Zt <- rho$Zt, "Matrix"))
    n <- nrow(X)
    p <- ncol(X)
    off <- rho$offset
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
        rho$RX <- if (is(rho$XtX <- crossprod(X), "dsparseMatrix"))
            Cholesky(rho$XtX - crossprod(rho$RZX)) else chol(rho$XtX)
    }
}

##' Install various objects, including Lambda, Lind, Zt and Ut in
##' environment rho based on the list bars of random effects terms and
##' model frame fr.
##'
##' @param bars a list of parsed random-effects terms
##' @param fr a model frame in which to evaluate these terms
##' @rho an environment that is modified by this function
##' @return NULL - the side effect of the function is to modify rho
makeZt <- function(bars, fr, rho) {
    stopifnot(is.list(bars), all(sapply(bars, is.language)),
              inherits(fr, "data.frame"), is.environment(rho))
    if (!length(bars))
        stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    blist <- lapply(bars,
                    function(x)
                {
                    ff <- eval(substitute(factor(fac), list(fac = x[[3]])), fr)
                    if (all(is.na(ff)))
                        stop("Invalid grouping factor specification,",
                             x[[3]])
                    nl <- length(levels(ff))
                    mm <- model.matrix(eval(substitute( ~ foo,
                                                       list(foo = x[[2]]))), fr)
                    nc <- ncol(mm)
                    nseq <- seq_len(nc)
                    sm <- as(ff, "sparseMatrix")
                    if (nc  > 1)
                        sm <- do.call(rBind, lapply(nseq, function(i) sm))
                    sm@x[] <- t(mm[])
                    ## When nc > 1 switch the order of the rows of sm
                    ## so the random effects for the same level of the
                    ## grouping factor are adjacent.
                    if (nc > 1)
                        sm <- sm[as.vector(matrix(seq_len(nc * nl),
                                                  nc = nl, byrow = TRUE)),]
                    list(ff = ff, sm = sm, nl = nl, cnms = colnames(mm))
                })
    nl <- sapply(blist, "[[", "nl")     # no. of levels per term
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
    nth <- as.integer((nc * (nc+1))/2)  # no. of parameters per term
    nb <- nc * nl                     # no. of random effects per term
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

    ## Need .f during transition of determinant(L) definition
    rho$.f <- if(package_version(packageDescription("Matrix")$Version) >=
                 "0.999375-31") 2 else 1
    NULL
}

lmer <-
    function(formula, data, REML = TRUE, sparseX = FALSE,
             control = list(), start = NULL, verbose = FALSE, doFit = TRUE,
             compDev = if (sparseX) FALSE else TRUE, subset, weights,
             na.action, offset, contrasts = NULL, ...)
{
    mf <- mc <- match.call()
    if (!is.null(list(...)$family)) {      # call glmer if family specified
        mc[[1]] <- as.name("glmer")
        eval(mc, parent.frame())
    }
    stopifnot(length(formula <- as.formula(formula)) == 3)
    if (missing(data)) data <- environment(formula)
                                        # environment for deviance evaluation
    rho <- new.env(parent = environment(formula))
    rho$REML <- as.logical(REML)[1]
                                        # evaluate and install the model frame
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula)
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
    } else as(model.matrix(fe.form, fr, contrasts), "dgeMatrix")
    rownames(X) <- NULL
    rho$X <- X
    rho$sparseX <- sparseX

    p <- ncol(X)
    stopifnot((rho$nmp <- rho$nobs - p) > 0)
    fixef <- numeric(p)
    names(fixef) <- colnames(X)
    rho$fixef <- fixef

    ## Check for method argument which is no longer used
    if (!is.null(method <- list(...)$method)) {
        msg <- paste("Argument", sQuote("method"),
                     "is deprecated.\nUse", sQuote("nAGQ"),
                     "to choose AGQ.  PQL is not available.")
        if (match.arg(method, c("Laplace", "AGQ")) == "Laplace") {
            warning(msg)
        } else stop(msg)
    }

    makeZt(expandSlash(findbars(formula[[3]])), fr, rho)
    rho$L <- Cholesky(tcrossprod(rho$Zt), LDL = FALSE, Imult = 1)
    derived_mats(rho)

    rho$mu <- numeric(rho$nobs)
    rho$gamma <- numeric(rho$nobs)
    rho$pwrss <- numeric(1)
    rho$ldL2 <- numeric(1)
    rho$ldRX2 <- numeric(1)
    rho$sqrtrwt <- sqrt(rho$weights)

    rho$compDev <- compDev
    rho$call <- mc
    
    sP <- function(x) {
### FIXME: weights are not yet incorporated (needed here?)
        if (compDev) {
            .Call("lmer_deviance", parent.env(environment()), x,
                  PACKAGE = "lme4a")
        } else {
            stopifnot(length(x) == length(theta))
            theta <<- as.numeric(x)
            Lambda@x[] <<- theta[Lind]           # update Lambda
            Ut <<- crossprod(Lambda, Zt)
            Matrix:::destructive_Chol_update(L, Ut, Imult = 1)
            cu <- solve(L, solve(L, crossprod(Lambda, Zty), sys = "P"),
                        sys = "L")
            RZX <<- solve(L, solve(L, crossprod(Lambda, ZtX), sys = "P"),
                          sys = "L")
            RX <<- if (sparseX) update(RX, XtX - crossprod(RZX)) else
                   chol(XtX - crossprod(RZX))
            fixef[] <<- if (sparseX) solve(RX, Xty - crossprod(RZX, cu))@x
            else solve(RX, solve(t(RX), Xty - crossprod(RZX, cu)))@x
            u[] <<- solve(L, solve(L, cu - RZX %*% fixef, sys = "Lt"),
                          sys = "Pt")@x
            gamma[] <<- (if(length(offset)) offset else 0) +
                (crossprod(Ut, u) + X %*% fixef)@x
            mu[] <<- gamma
            pwrss <<- sum(c(y - mu, u)^2) # penalized residual sum of squares
            ldL2[] <<- .f * determinant(L)$mod
            ldRX2[] <<- 2 * determinant(RX)$mod
            if (!REML) return(ldL2 + nobs * (1 + log(2 * pi * pwrss/nobs)))
            ldL2 + ldRX2 + nmp * (1 + log(2 * pi * pwrss/nmp))
        }
    }
    gP <- function() theta
    gB <- function() cbind(lower = lower,
                           upper = rep.int(Inf, length(theta)))
    environment(sP) <- environment(gP) <- environment(gB) <- rho
    me <- new("lmerenv", setPars = sP, getPars = gP, getBounds = gB)
    if (doFit) {                        # perform the optimization
### FIXME: Allow for selecting an optimizer.  Use optimize for scalar
### problems
        if (verbose) control$trace <- 1
        nlminb(is.finite(rho$lower), me@setPars, lower = rho$lower,
               control = control)
    }
    me
}

glmer <-
function(formula, data, family = gaussian, sparseX = FALSE,
         compDev = if (sparseX) FALSE else TRUE,
         control = list(), start = NULL, verbose = FALSE, doFit = TRUE,
         subset, weights, na.action, offset, contrasts = NULL, nAGQ = 1,
         mustart, etastart, ...)
{
    mf <- mc <- match.call()
    if (missing(family)) {
        mc[[1]] <- as.name("lmer")
        eval(mc, parent.frame())
    }

    ## Check for invalid specifications
    if (!is.null(method <- list(...)$method)) {
        msg <- paste("Argument", sQuote("method"),
                     "is deprecated.\nUse", sQuote("nAGQ"),
                     "to choose AGQ.  PQL is not available.")
        if (match.arg(method, c("Laplace", "AGQ")) == "Laplace") {
            warning(msg)
        } else stop(msg)
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
    fixef <- numeric(p)
    names(fixef) <- colnames(X)
    rho$fixef <- fixef
    rho$start <- numeric(p);            # needed for family$initialize

    ## evaluate, check and initialize family
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if(is.function(family)) family <- family()
    rho$family <- family
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

    ## enforce modes
    rho$y <- unname(as.double(rho$y))   # must be done after initialize
    rho$pwrss <- numeric(1)
    rho$ldL2 <- numeric(1)
    rho$ldRX2 <- numeric(1)
    rho$sqrtrwt <- numeric(n)

    makeZt(expandSlash(findbars(formula[[3]])), fr, rho)
    rho$L <- Cholesky(tcrossprod(rho$Zt), LDL = FALSE, Imult = 1)
    derived_mats(rho)
    rho$compDev <- compDev
    rho
}

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
                   verbose = FALSE, nAGQ = 1, doFit = TRUE, subset,
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
    ft <- famType(family)
    rho$dims[names(ft)] <- ft
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
    rho$start[] <- rho$fixef <- qr.coef(qrX, unlist(lapply(pnames, get, envir = rho$nlenv)))
    rho$RX <- qr.R(qrX)
    eb <- evalbars(formula, fr, contrasts, TRUE) # flist, trms, nest
    rho$dims["nest"] <- eb$nest
    rho$flist <- eb$flist
    lmerFactorList(eb$trms, rho)

    q <- length(rho$u)
    rho$u0 <- numeric(q)
                                        # evaluate the control argument
    control <- do.call(lmerControl, as.list(control))
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
    rho$beta0[] <- rho$fixef
    rho$u0[] <- rho$u
    if (!doFit) return(rho)

#    merFinalize(rho)
}

setMethod("fixef", "merenv", function(object, ...) env(object)$fixef)

setMethod("ranef", "merenv", function(object, ...)
          {
              with(env(object), (Lambda %*% u)@x)
          })

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
          nc <- sapply(cnms <- rho$cnms, length)
          nb <- nc * (nl <- unname(sapply(levs, length))[asgn])
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

devcomp <- function(x, ...) {
    stopifnot(is(x, "lmerenv"))
    with(env(x),
         list(cmp = c(ldL2 = ldL2, ldRX2 = ldRX2, pwrss = pwrss,
              deviance = ldL2 + nobs * (1 + log(2 * pi * pwrss/nobs)),
              REML = ldL2 + ldRX2 + nmp * (1 + log(2 * pi * pwrss/nmp))),
              dims = c(n = nobs, p = length(fixef), nmp = nmp, q = nrow(Zt))))
}

setMethod("sigma", signature(object = "lmerenv"),
          function (object, ...) {
              dc <- devcomp(object)
              unname(sqrt(dc$cmp["pwrss"]/
                          dc$dims[if (env(object)$REML) "nmp" else "n"]))
          })

## Obtain the ML fit of the model parameters
MLfit <- function(x) {
    stopifnot(is(x, "lmerenv"))
    if (!env(x)$REML) return(x)
    update(x, REML = FALSE)
}
    
setMethod("anova", signature(object = "lmerenv"),
	  function(object, ...)
      {
	  mCall <- match.call(expand.dots = TRUE)
	  dots <- list(...)
	  modp <- as.logical(sapply(dots, is,
                                    "lmerenv")) | as.logical(sapply(dots,
                                                                    is, "lm"))
	  if (any(modp)) {		# multiple models - form table
	      opts <- dots[!modp]
	      mods <- c(list(object), dots[modp])
	      names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)],
				    as.character)
              mods <- lapply(mods, function(x) ## get ML estimates if necessary
                         {
                             if (!env(x)$REML) return(x)
                             update(x, REML = FALSE)
                         })
	      mods <- mods[order(sapply(lapply(mods, logLik), attr, "df"))]
	      llks <- lapply(mods, logLik)
	      Df <- sapply(llks, attr, "df")
	      calls <- lapply(mods, function(mod) env(mod)$call)
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
				AIC = sapply(llks, AIC),
				BIC = sapply(llks, BIC),
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
              p <- dc$dims["p"]
              ss <- fixef(object)
              stop("assign attribute not currently available")
              asgn <- attr(object@X, "assign")
              terms <- terms(object)
              nmeffects <- attr(terms, "term.labels")
              if ("(Intercept)" %in% names(ss))
                  nmeffects <- c("(Intercept)", nmeffects)
              ss <- unlist(lapply(split(ss, asgn), sum))
              df <- unlist(lapply(split(asgn,  asgn), length))
              ## dfr <- unlist(lapply(split(dfr, asgn), function(x) x[1]))
              ms <- ss/df
              f <- ms/(sigma(object)^2)
              ## P <- pf(f, df, dfr, lower.tail = FALSE)
              ## table <- data.frame(df, ss, ms, dfr, f, P)
              table <- data.frame(df, ss, ms, f)
              dimnames(table) <-
                  list(nmeffects,
                                        #			c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
                       c("Df", "Sum Sq", "Mean Sq", "F value"))
              if ("(Intercept)" %in% nmeffects)
                  table <- table[-match("(Intercept)", nmeffects), ]
              attr(table, "heading") <- "Analysis of Variance Table"
              class(table) <- c("anova", "data.frame")
              table
	  }
      })


##' Extract the conditional variance-covariance matrix of the fixed
##' effects
setMethod("vcov", signature(object = "lmerenv"),
	  function(object, ...)
      {
          en <- env(object)
          rr <- as(sigma(object)^2 * chol2inv(en$RX), "dpoMatrix")
          nms <- colnames(en$X)
          dimnames(rr) <- list(nms, nms)
          rr@factors$correlation <- as(rr, "corMatrix")
          rr
      })

##' Create the VarCorr object of variances and covariances
setMethod("VarCorr", signature(x = "lmerenv"),
	  function(x, ...)
      {
          sc <- sigma(x)
          rho <- env(x)
          nc <- sapply(cnms <- rho$cnms, length)
          ncseq <- seq_along(nc)
          thl <- split(rho$theta, rep.int(ncseq, (nc * (nc + 1))/2))
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
          fl <- rho$flist
          names(ans) <- names(fl)[attr(fl, "assign")]
          attr(ans, "sc") <- sc
          ans
      })

## This is modeled a bit after  print.summary.lm :
printMerenv <- function(x, digits = max(3, getOption("digits") - 3),
                        correlation = TRUE, symbolic.cor = FALSE,
                        signif.stars = getOption("show.signif.stars"), ...)
{
    so <- summary(x)
    cat(so$methTitle, "\n")
    if (!is.null(so$call$formula))
        cat("Formula:", deparse(so$call$formula),"\n")
    if (!is.null(so$call$data))
        cat("   Data:", deparse(so$call$data), "\n")
    if (!is.null(so$call$subset))
        cat(" Subset:",
            deparse(asOneSidedFormula(x$call$subset)[[2]]),"\n")
    print(so$AICtab, digits = digits)
    cat("\nRandom effects:\n")
    print(so$REmat, quote = FALSE, digits = digits, ...)

    ngrps <- so$ngrps
    cat(sprintf("Number of obs: %d, groups: ", so$devcomp$dims["n"]))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
    if (nrow(so$coefficients) > 0) {
	cat("\nFixed effects:\n")
	printCoefmat(so$coefficients, zap.ind = 3, #, tst.ind = 4
		     digits = digits, signif.stars = signif.stars)
	if(correlation) {
	    corF <- so$vcov@factors$correlation
	    if (!is.null(corF)) {
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
}

setMethod("summary", signature(object = "lmerenv"),
	  function(object, ...)
      {
          dc <- devcomp(object)
          rho <- env(object)
          REML <- rho$REML
          fcoef <- fixef(object)
          vcov <- vcov(object)
          corF <- vcov@factors$correlation
          coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF@sd)
          llik <- logLik(object, REML)
          mName <- paste("Linear mixed model fit by",
                         if (REML) "REML" else "maximum likelihood")
          AICframe <- if (REML) dc$cmp["REML"] else
          data.frame(AIC = AIC(llik), BIC = BIC(llik),
                     logLik = c(llik),
                     deviance = unname(dc$cmp["deviance"]), row.names = "")
          varcor <- VarCorr(object)
          REmat <- formatVC(varcor)
          if (nrow(coefs) > 0)
              coefs <- cbind(coefs, "t value" = coefs[,1]/coefs[,2])
          ans <- list(methTitle = mName,
                      devcomp = dc,
                      logLik = llik,
                      ngrps = sapply(rho$flist, function(x) length(levels(x))),
                      sigma = sigma(object),
                      coefficients = coefs,
                      vcov = vcov,
                      REmat = REmat,
                      AICtab= AICframe,
                      call = rho$call
                      )
          class(ans) <- "summary.lmer" # use S3 class for now
          ans
      })## summary()

setMethod("model.frame", signature(formula = "lmerenv"),
	  function(formula, ...) env(formula)$frame)

setMethod("model.matrix", signature(object = "lmerenv"),
	  function(object, ...) env(object)$X)

setMethod("terms", signature(x = "lmerenv"),
	  function(x, ...) attr(env(x)$frame, "terms"))

setMethod("deviance", signature(object="lmerenv"),
	  function(object, REML = NULL, ...)
      {
          if (missing(REML) || is.null(REML) || is.na(REML[1]))
              REML <- env(object)$REML
          devcomp(object)$cmp[if(REML) "REML" else "deviance"]
      })

##' Extract the log-likelihood or restricted log-likelihood
setMethod("logLik", signature(object="lmerenv"),
	  function(object, REML = NULL, ...)
      {
          rho <- env(object)
          if (is.null(REML) || is.na(REML[1]))
              REML <- rho$REML
          val <- -deviance(object, REML = REML)/2
          attr(val, "nall") <- attr(val, "nobs") <- length(rho$y)
          attr(val, "df") <- length(rho$fixef) +
              length(rho$theta) + 1L
          attr(val, "REML") <-  as.logical(REML)
          class(val) <- "logLik"
          val
      })

setMethod("update", signature(object = "merenv"),
	  function(object, formula., ..., evaluate = TRUE)
      {
	  call <- env(object)$call
	  if (is.null(call))
	      stop("env(object) should contain a call object")
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
      })

setMethod("print", "lmerenv", printMerenv)
setMethod("show", "lmerenv", function(object) printMerenv(object))

setMethod("coef", signature(object = "merenvtrms"),
	  function(object, ...)
      {
          if (length(list(...)))
              warning(paste('arguments named "',
                            paste(names(list(...)), collapse = ", "),
                                  '" ignored', sep = ''))
          fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
          ref <- ranef(object)
          val <- lapply(ref, function(x) fef[rep(1, nrow(x)),,drop = FALSE])
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
       })

