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
    nl <- unlist(lapply(blist, "[[", "nl"))# no. of levels per term

    ## order terms stably by decreasing number of levels in the factor
    if (any(diff(nl)) > 0) {
        ord <- rev(order(nl))
        blist <- blist[ord]
        nl <- nl[ord]
    }
    Zt <- do.call(rBind, lapply(blist, "[[", "sm"))
    q <- nrow(Zt)

    ## Create and install Lambda, Lind, etc.  This must be done after
    ## any potential reordering of the terms.
    cnms <- lapply(blist, "[[", "cnms")
    nc <- sapply(cnms, length)          # no. of columns per term
    nth <- as.integer((nc * (nc+1))/2)  # no. of parameters per term
    nb <- nc * nl                       # no. of random effects per term
    stopifnot(sum(nb) == q)
    thet <- numeric(sum(nth))

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
    ll <- list(Class = "reTrms", Zt = Zt, u = numeric(q), ##-nL nLevs = nl,
               theta = thet, Lind = as.integer(Lambda@x))
    ## lower bounds on theta elements are 0 if on diagonal, else -Inf
    ll$lower <- -Inf * (thet + 1)
    ll$lower[unique(diag(Lambda))] <- 0
    ll$theta[] <- is.finite(ll$lower)   # initial values of theta are 0 off-diagonal, 1 on
    Lambda@x[] <- ll$theta[ll$Lind]     # initialize elements of Lambda
    ll$Lambda <- Lambda

    Ut <- if (s > 1) { ## Ut is the sum of vertical sections of Zt
        N <- ncol(Zt)
        Reduce("+",
               lapply(split(seq_len(N),
                            rep.int(seq_len(s),
                                    rep.int(N %/% s, s))),
                      function(cols) Zt[, cols]))
    }
    else crossprod(Lambda, Zt)

    ll$L <- Cholesky(tcrossprod(Ut), LDL = FALSE, Imult = 1)
                                        # massage the factor list
    fl <- lapply(blist, "[[", "ff")
                                        # check for repeated factors
    fnms <- names(fl)
    if (length(fnms) > length(ufn <- unique(fnms))) {
        fl <- fl[match(ufn, fnms)]
        asgn <- match(fnms, ufn)
    } else asgn <- seq_along(fl)
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
##'
##' @return an object that inherits from feModule
mkFeModule <- function(form, fr, contrasts, reMod, sparseX)
{
    ## fixed-effects model matrix X - remove random parts from formula:
    form[[3]] <- if(is.null(nb <- nobars(form[[3]]))) 1 else nb
    X <- model.Matrix(form, fr, contrasts, sparse = sparseX,
		      row.names = FALSE)
    ZtX <- reMod@Zt %*% X
    new(if(sparseX) "spFeMod" else "deFeMod",
	RZX = ZtX,
	X = X,
	coef = numeric(ncol(X)),
	fac = {
	    if (sparseX)    # crossprod(ZtX) may be more dense than X.X
		Cholesky(crossprod(X) + crossprod(ZtX), LDL = FALSE)
	    else
		chol(crossprod(X)) })
}

##' <description>
##' Generic function to update the dimensions component of the devcomp
##' slot in a merMod object.  The devcomp slot is passed from module
##' to module accumulating information along the way.
##' <details>
##' @title Update the dimensions in the devcomp slot
##' @param x object on which the update is based
##' @param dcmp current value of the devcomp slot
##' @return updated devcomp slot
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
              dcmp$dims["p"] <- length(x@coef)
              dcmp
          })
setMethod("updateDcmp", signature(x = "spFeMod", dcmp = "list"),
          function(x, dcmp) {
              dcmp$dims["spFe"] <- 1L
              dcmp$dims["p"] <- length(x@coef)
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
setMethod("updateDcmp", signature(x = "glmRespMod", dcmp = "list"),
          function(x, dcmp) {
              dcmp <- .respBase(x, dcmp)
              dcmp$dims["GLMM"] <- 1L
              dcmp$dims["useSc"] <- 0L
              dcmp
          })
setMethod("updateDcmp", signature(x = "nlsRespMod", dcmp = "list"),
          function(x, dcmp) {
              dcmp <- .respBase(x, dcmp)
              dcmp$dims["NLMM"] <- 1L
              dcmp
          })
setMethod("updateDcmp", signature(x = "nglmRespMod", dcmp = "list"),
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

lmer <- function(formula, data, REML = TRUE, sparseX = FALSE,
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
    dcmp <- updateDcmp(reTrms, .dcmp())
                                        # fixed-effects module
    feMod <- mkFeModule(formula, fr, contrasts, reTrms, sparseX)
    rr <- mkRespMod(fr)
    class(rr) <- "lmerResp"
    rr@REML <- ifelse(REML, ncol(feMod@X), 0L)
    ans <- new("merMod",
               call = mc,
               devcomp = updateDcmp(rr, updateDcmp(feMod, dcmp)),
               frame = fr,
	       re = reTrms, fe = feMod, resp = rr)
    if (doFit) {                        # optimize estimates
        if (verbose) control$iprint <- 2L
        devfun <- mkdevfun(ans, compDev = compDev)
        opt <- bobyqa(ans@re@theta, devfun, ans@re@lower, control = control)
        ans <- updateMod(ans, opt$par, opt$fval)
    }
    ans
}## { lmer }

##' <description>
##' Update the components of an merMod object from the results of an
##' optimization.  To maintain proper semantics for an R function the
##' object passed to the deviance evaluation is treated as read-only.
##' Once the deviance is optimized the model object is updated with
##' this function.
##' <details>
##' @title Update the components of an merMod object.
##' @param mod an merMod object
##' @param pars values of the parameters at which to update
##' @param fval deviance evaluation at the estimated parameters
##'   \code{pars}.  Note that this deviance evaluation is a numeric
##'  vector of length 1 with several attributes.  If the optimizer
##' used does not return this evaluation with its attributes you will
##' need to evaluate the deviance function before calling
##' \code{updateMod}.  (We could use the parameters and the deviance
##' function as arguments but that would entail an extra deviance
##' function evaluation, which can be expensive for exotic models.)
##' @return the updated merMod object
updateMod <- function(mod, pars, fval) {
    stopifnot(is(mod, "merMod"),
              is.numeric(pars),
              is.numeric(fval),
              is.numeric(u <- attr(fval, "u")),
              is.numeric(beta <- attr(fval, "beta")),
              is.numeric(ldL2 <- attr(fval, "ldL2")),
              length(ldL2) == 1L,
              is.numeric(wrss <- attr(fval, "wrss")),
              length(wrss) == 1L,
              is.numeric(ussq <- attr(fval, "ussq")),
              length(ussq) == 1L)
    re <- mod@re
    fe <- mod@fe
    resp <- mod@resp
    devcomp <- mod@devcomp
    lower <- re@lower
    stopifnot(length(pars) >= length(lower),
              length(beta) == ncol(fe@X),
              length(u) == nrow(re@Zt))
    pars <- pars[seq_along(lower)]

    re@Lambda@x[] <- pars[re@Lind]      # update Lambda
    re@theta <- pars
    re@u <- u
    ## FIXME: Create the reModule as inheriting from respModule with Z
    ## %*% Lambda as the model matrix X and u as coef.  Then the only
    ## thing that needs to be changed is to allow the multiple of the
    ## identity in the update of fac.
    resp <- updateWts(updateMu(resp, as.vector(crossprod(re@Zt, re@Lambda %*% u)
                                               + fe@X %*% beta)))
    
    fe@coef <- beta
    fe <- reweightPred(fe, resp@sqrtXwt, resp@wtres)
    ## reweight the re module.  This should eventually be a method.
    WtMat <- Diagonal(x=as.vector(resp@sqrtXwt))
    ## FIXME: the previous calculation will not work for an nlsRespMod
    Ut <- crossprod(re@Lambda, re@Zt %*% WtMat)
    re@L <- update(re@L, Ut, 1)
    V <- WtMat %*% fe@X
    UtV <- Ut %*% V
    fe@RZX <- solve(re@L, solve(re@L, UtV, sys = "P"), sys = "L")
    fe@fac <- if (is(V, "sparseMatrix")) {
        Cholesky(crossprod(V) - crossprod(fe@RZX), LDL = FALSE) # update instead?
    } else {
        chol(crossprod(V) - crossprod(fe@RZX))
    }
    dc <- mod@devcomp$cmp
    dc["ldRX2"] <- 2 * determinant(fe@fac, log = TRUE)$modulus
    dc[ifelse(is(resp, "lmerResp") && resp@REML, "REML", "dev")] <- as.vector(fval)
    dc["ldL2"] <- ldL2
    dc["wrss"] <- wrss
    dc["ussq"] <- ussq
    dc["pwrss"] <- pwrss <- wrss + ussq
    dc["sigmaML"] <- sqrt(pwrss/length(resp@y))
    dc["sigmaREML"] <- sqrt(pwrss/(length(resp@y)-length(beta)))
    mod@devcomp$cmp <- dc
    mod@re <- re
    mod@fe <- fe
    mod@resp <- resp
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
mkdevfun <- function(mod, nAGQ = 1L, u0 = numeric(length(mod@re@u)), verbose = 0L, compDev = TRUE) {
    stopifnot(is(mod, "merMod"))
    resp <- mod@resp
    beta0 <- numeric(length(mod@fe@coef))
    if (is(resp, "lmerResp")) {
        if (compDev) return(function(th) .Call("merDeviance", mod, th, beta0, u0, 0L, 3L, PACKAGE="lme4a"))
##        if (compDev) return(function(th) .Call(merDeviance, mod, th, beta0, u0, 0L, 3L)
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
            wrss <- sum(((resp@y - mu)*resp@sqrtrwt)^2)
            ussq <- sum(u@x^2)
            pwrss <- wrss + ussq
            ldL2 <- as.vector(2 * determinant(L)$mod)
            if (resp@REML) {
                nmp <- length(mu) - length(beta)
                ldRX2 <- as.vector(2 * determinant(RX)$mod)
                ans <- ldL2 + ldRX2 + nmp * (1 + log(2 * pi * pwrss/nmp))
            } else {
                n <- length(mu)
                ans <- ldL2	     + n   * (1 + log(2 * pi * pwrss/n  ))
            }
            structure(unname(ans), beta = beta, u = u@x, ldL2 = ldL2,
                      wrss = wrss, ussq = ussq)
        })
    }
    if (nAGQ == 0L)
#        return(function(pars) .Call(merDeviance, mod, pars, beta0, u0, verbose, 3L))
        return(function(pars) .Call("merDeviance", mod, pars, beta0, u0, verbose, 3L, PACKAGE="lme4a"))
    thpars <- seq_along(mod@re@theta)
#    function(pars) .Call(merDeviance, mod, pars[thpars], pars[-thpars], u0, verbose, 2L)
    function(pars) .Call("merDeviance", mod, pars[thpars], pars[-thpars], u0, verbose, 2L, PACKAGE="lme4a")
}

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
	  as.vector(X %*% object@fe @ coef) +  # fixed-effect contribution
	      sigma(object) * (## random-effects contribution:
			       if(use.u) {
				   object@re @ u
			       } else {
				   U <- crossprod(object@re@Zt, object@re@Lambda)
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
        U <- crossprod(x@re @ Zt, x@re @ Lambda)
        q <- ncol(U)
    }
    Zt <- x@re @ Zt
    X.beta <- as.vector(X %*% x@fe @ coef) # fixed-effect contribution
    sigm.x <- sigma(x)

    ## Here, and below ("optimize"/"bobyqa") using the "logic" of lmer2() itself:
## lmer..Update <- if(is(x, "lmerSp")) lmerSpUpdate else lmerDeUpdate
#    devfun <- mkdevfun(x)
##    oneD <- length(x@re@theta) < 2
    theta0 <- x@re@theta
    ## just for the "boot" result -- TODOmaybe drop
    mle <- list(beta = x@fe @ coef, theta = theta0, sigma = sigm.x)

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
        opt <- bobyqa(theta0, mkdevfun(x), x@re@lower, control = control)
        xx <- updateMod(x, opt$par, opt$fval)
            ## FIXME: also here, prefer \hat\sigma^2 == 0 (exactly)
##        }
        foo <- tryCatch(FUN(xx), error = function(e)e)
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
    opt <- bobyqa(ans@re@theta, mkdevfun(ans, 0L, verbose = verbose),
                  ans@re@lower, control = control)
    if (nAGQ > 0L) {
        thpars <- seq_along(ans@re@theta)
        bb <- attr(opt$fval, "beta")
        opt <- bobyqa(c(opt$par, bb),
                      mkdevfun(ans, nAGQ, attr(opt$fval, "u"), verbose),
                      lower = c(ans@re@lower, rep.int(-Inf, length(bb))),
                      control = control)
    }
    updateMod(ans, opt$par, opt$fval)
}

glmer <- function(formula, data, family = gaussian, sparseX = FALSE,
                  control = list(), start = NULL, verbose = 0L, nAGQ = 1L,
                  doFit = TRUE, subset, weights, na.action, offset,
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

    feMod <- mkFeModule(formula, fr, contrasts, reTrms, sparseX)
    respMod <- mkRespMod(fr, family)
    ans <- new("merMod",
               call = mc,
               devcomp = updateDcmp(respMod, updateDcmp(feMod, dcmp)),
               frame = fr,
	       re = reTrms, fe = feMod, resp = respMod)
    if (!doFit) return(ans)
    ans <- PIRLSest(ans, verbose, control, nAGQ)

}## {glmer}

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
nlmer <- function(formula, data, family = gaussian, start = NULL,
                  verbose = 0L, nAGQ = 1L, doFit = TRUE,
                  subset, weights, na.action, mustart, etastart,
                  sparseX = FALSE, contrasts = NULL, control = list(), ...)
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
    feMod <- mkFeModule(fe.form, frE, contrasts, reTrms, sparseX = sparseX)
                                        # should this check be in mkFeModule?
    p <- length(feMod@coef)
    if ((qrX <- qr(feMod@X))$rank < p)
        stop(gettextf("rank of X = %d < ncol(X) = %d", qrX$rank, p))
    feMod@coef <- qr.coef(qrX, unlist(lapply(pnames, get, envir = nlenv)))
    respMod <- mkRespMod(fr, nlenv = nlenv, nlmod = nlmod)
    ans <- new("merMod",
               call = mc,
               devcomp = updateDcmp(respMod, updateDcmp(feMod, dcmp)),
               frame = fr, re = reTrms, fe = feMod, resp = respMod)
    if (!doFit) return(ans)
    PIRLSest(ans, verbose, control, nAGQ)
}

## Methods for the merMod class
setMethod("fixef",  "merMod", function(object, ...)
          structure(object@fe@coef, names = dimnames(object@fe@X)[[2]]))

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
#                  vv <- .Call(reTrmsCondVar, re, sigma(object))
                  vv <- .Call("reTrmsCondVar", re, sigma(object), PACKAGE="lme4a")
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

##' update()  for all kind  "merMod", "merenv", "mer", .. :
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

setMethod("update", signature(object = "merMod"), updateMer)

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
    tab <- so$AICtab
    if (length(tab) == 1 && names(tab) == "REML")
	cat("REML criterion at convergence:", round(tab, 4), "\n")
    else print(round(so$AICtab, 4))
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
		cat(sprintf(paste("\nCorrelation matrix not shown by default, as p = %d > 20.",
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

setMethod("print", "merMod", printMerenv)
setMethod("show",  "merMod", function(object) printMerenv(object))
setMethod("fitted", "merMod", function(object, ...) {object <- object@resp; callGeneric(...)})
#setMethod("fitted", "merMod", function(object,...) object@resp@mu)
setMethod("residuals", "merMod",
          function(object, type = c("deviance", "pearson",
                           "working", "response", "partial"), ...)
      {
          object <- object@resp
          callGeneric(...)
      })

setMethod("print", "summary.mer", printMerenv)
setMethod("show",  "summary.mer", function(object) printMerenv(object))


## coef() method for all kinds of "mer", "*merMod", ... objects
## ------  should work with fixef() + ranef()  alone
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

setMethod("coef", signature(object = "merMod"), coefMer)

## FIXME: Do we really need a separate devcomp extractor?  I suppose it can't hurt.
setMethod("devcomp", "merMod", function(x, ...) x@devcomp)

setMethod("getL", "reModule", function(x) x@L)
setMethod("getL", "merMod", function(x) x@re@L)

setMethod("isREML", "merMod", function(x) as.logical(x@devcomp$dims["REML"]))

setMethod("refitML", "merMod",
          function (x) {
              if (!isREML(x)) return(x)
              x@resp@REML <- 0L
              x@devcomp$dims["REML"] <- 0L
              opt <- bobyqa(x@re@theta, mkdevfun(x), x@re@lower)
              updateMod(x, opt$par, opt$fval)
          })

setMethod("getCall", "merMod",	function(x) x@call)

## A more effective method for merMod objects is defined in lmer2.R

##' <description>
##'
##' <details>
##' @title anova() for both  "lmerenv" and "lmer" fitted models
##' @param object an "lmerenv", "lmer" or "lmerMod" - fitted model
##' @param ...  further such objects
##' @return an "anova" data frame; the traditional (S3) result of anova()
anovaLmer <- function(object, ...) {
    mCall <- match.call(expand.dots = TRUE)
    dots <- list(...)
    .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))
    modp <- {
#	as.logical(.sapply(dots, is, "lmerenv")) |
	as.logical(.sapply(dots, is, "merMod")) |
#	as.logical(.sapply(dots, is, "lmer")) |
	as.logical(.sapply(dots, is, "lm")) }
    if (any(modp)) {			# multiple models - form table
	opts <- dots[!modp]
	mods <- c(list(object), dots[modp])
	## model names
	mNms <- .sapply(as.list(mCall)[c(FALSE, TRUE, modp)], deparse)
	names(mods) <- sub("@env$", '', mNms) # <- hack
	mods <- lapply(mods, refitML)

        devs <- sapply(mods, deviance)
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
                          deviance = -2*llk,
			  Chisq = chisq,
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
	asgn <- object@fe @ X @ assign
	stopifnot(length(asgn) == p,
                  is(object@fe, "deFeMod"), # haven't worked out sparse version
                  is(object@re, "reTrms"))  # things are really weird with no re terms
        ss <- ((object@fe@fac %*% object@fe@coef)@x)^2
        names(ss) <- colnames(object@fe@X)
        terms <- terms(object)
        nmeffects <- setdiff(attr(terms, "term.labels"), names(object@re@flist))
        if ("(Intercept)" %in% names(ss))
	    nmeffects <- c("(Intercept)", nmeffects)
	ss <- unlist(lapply(split(ss, asgn), sum))
        stopifnot(length(ss) == length(nmeffects))
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
}## {anovaLmer}

setMethod("anova", signature(object = "merMod"), anovaLmer)

##' <description>
##'
##' <details>
##' @title vcov(): Extract conditional covariance matrix of fixed effects
##' @param sigma = sigma(object)
##' @param correlation
##' @param ...
mkVcov <- function(sigma, RX, nmsX, correlation = TRUE, ...) {
    V <- sigma^2 * chol2inv(RX)
    if(is.null(rr <- tryCatch(as(V, "dpoMatrix"),
                              error = function(e) NULL)))
        stop("Computed variance-covariance matrix is not positive definite")
    dimnames(rr) <- list(nmsX, nmsX)
    if(correlation)
        rr@factors$correlation <- as(rr, "corMatrix")
    rr
}

setMethod("vcov", signature(object = "merMod"),
	  function(object, correlation = TRUE, sigm = sigma(object), ...)
	  mkVcov(sigm, RX = object@fe@fac, nmsX = colnames(object@fe@X),
		 correlation=correlation, ...))

setMethod("vcov", "summary.mer",
	  function(object, correlation = TRUE, ...)
      {
	  if(!is.null(object$vcov))
	      vcov
	  else if(!is.null(FE <- object$fe))
	      mkVcov(object$sigma, RX = FE@fac, nmsX = colnames(FE@X),
		     correlation=correlation, ...)
	  else stop("Both 'vcov' and 'fe' components are missing.  You need\n",
		    "at least one TRUE in summary(..,  varcov = *, keep.X = *)")
      })


### "FIXME": instead of 'Lambda', it is sufficient
##  ------   to just keep a list of the small lower triangular blocks
## 'Li' from above;  then, Lambda = bdiag({<Li>}), i.e.  do.call(bdiag, blocksLambda(..))
blocksLambda <- function(re) {
    if (!is(re, "reTrms")) stop("only works for \"reTrms\"")
    cnms <- re@cnms
    nc <- unlist(lapply(cnms, length),recursive=FALSE) # no. of columns per term
    ncseq <- seq_along(nc)
    thl <- split(re@theta, rep.int(ncseq, (nc * (nc + 1))/2))
    fl <- re@flist
    structure(lapply(ncseq, function(i)
		 {
		     ## Li := \Lambda_i, the i-th block diagonal of \Lambda(\theta)
		     Li <- diag(nrow = nc[i])
		     Li[lower.tri(Li, diag = TRUE)] <- thl[[i]]
		     rownames(Li) <- cnms[[i]]
		     Li
		 }),
	      names = names(fl)[attr(fl, "assign")],
	      ncols = nc,
              nLevs = sapply(fl, nlevels)
	      ##-nL nLevs = re@nLevs
              )
}
if(FALSE) { ## Now we can easily "build" Lambda (inside the 're'):
    re <- mod@re
    bL <- lme4a:::blocksLambda(re)
    print( sapply(bL, rcond) ) ## reciprocal condition numbers
    identical(re@Lambda, as(Matrix:::.bdiag(rep.int(bL, attr(bL,"nLevs"))),
			    "CsparseMatrix"))
}
setMethod("rcond", signature(x = "reTrms", norm = "character"),
          function(x, norm, ...)
          sapply(blocksLambda(x), rcond, norm=norm))


## Keep this separate, as it encapsulates the computation
## w/o explicit use of S4 modules
mkVarCorr <- function(sc, cnms, nc, theta, nms) {
    ncseq <- seq_along(nc)
    thl <- split(theta, rep.int(ncseq, (nc * (nc + 1))/2))
    ans <- lapply(ncseq, function(i)
              {
		  ## Li := \Lambda_i, the i-th block diagonal of \Lambda(\theta)
		  Li <- diag(nrow = nc[i])
		  Li[lower.tri(Li, diag = TRUE)] <- thl[[i]]
		  rownames(Li) <- cnms[[i]]
		  ## val := \Sigma_i = \sigma^2 \Lambda_i \Lambda_i', the
		  val <- tcrossprod(sc * Li) # variance-covariance
                  stddev <- sqrt(diag(val))
                  correl <- t(val / stddev)/stddev
                  diag(correl) <- 1
                  attr(val, "stddev") <- stddev
                  attr(val, "correlation") <- correl
                  val
              })
    if(is.character(nms)) names(ans) <- nms
    attr(ans, "sc") <- sc
    ans
}

setMethod("VarCorr", signature(x = "merMod"),
          function(x, sigma, rdig)# <- 3 args from nlme
      {
	  if (!is(re <- x@re, "reTrms"))
	      stop("VarCorr methods require reTrms, not just reModule")
	  cnms <- re@cnms
	  if(missing(sigma)) # "bug": fails via default 'sigma=sigma(x)'
	      sigma <- sigma(x)
	  nc <- sapply(cnms, length)	  # no. of columns per term
	  mkVarCorr(sigma, cnms=cnms, nc=nc, theta = re@theta,
		    nms = {fl <- re@flist; names(fl)[attr(fl, "assign")]})
      })

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
unscaledVar <- function(object, RX = object@fe@fac)
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
			   function(x) {
			       x <- as(x, "matrix")
			       cc <- format(round(x, 3), nsmall = 3)
			       cc[!lower.tri(cc)] <- ""
			       nr <- dim(cc)[1]
			       if (nr >= maxlen) return(cc)
			       cbind(cc, matrix("", nr, maxlen-nr))
			   }))[, -maxlen, drop = FALSE]
        if (nrow(corr) < nrow(reMat))
            corr <- rbind(corr, matrix("", nr = nrow(reMat) - nrow(corr), nc = ncol(corr)))
	colnames(corr) <- rep.int("", ncol(corr))
        colnames(corr)[1] <- "Corr"
	cbind(reMat, corr)
    } else reMat
}

##' <description>
##'
##' @title Summary Method for *mer() fits, i.e., "merMod" objects
##' @param object
##' @param varcov logical indicating if vcov(.) should be computed and stored.
##' @param keep.X logical indicating if the 'fe' component of object should be stored;
##'   the default is true when 'varcov' is false, as we then need fe for vcov()
##' @param ...
##' @return S3 class "summary.mer", basically a list .....
setMethod("summary", "merMod",
          function(object, varcov = FALSE, keep.X = !varcov, ...)
      {
          resp <- object@resp
          devC <- object@devcomp
          dd <- devC$dims
          cmp <- devC$cmp
          useSc <- as.logical(dd["useSc"])
          sig <- sigma(object)
          REML <- isREML(object)

          fam <- NULL
          if(is(resp, "glmRespMod")) fam <- resp@family
          coefs <- cbind("Estimate" = fixef(object),
                         "Std. Error" = sig * sqrt(unscaledVar(RX = object@fe@fac)))
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
	  ## FIXME: You can't count on object@re@flist,
	  ##        nor compute VarCorr() unless is(re, "reTrms"):
	  varcor <- VarCorr(object)
                                        # use S3 class for now
          structure(list(methTitle = mName,
                         devcomp = devC, isLmer = is(resp, "lmerResp"), useScale = useSc,
                         logLik = llik, family = fam,
			 ngrps = sapply(object@re@flist, function(x) length(levels(x))),
                         coefficients = coefs,
                         sigma = sig,
                         vcov = if(varcov) vcov(object, correlation=TRUE, sigm=sig),
                         fe = if(keep.X) object@fe,
                         varcor = varcor, # and use formatVC(.) for printing.
                         AICtab= AICstats,
                         call = object@call
                         ), class = "summary.mer")
      })

setMethod("summary", "summary.mer",
	  function(object, varcov = FALSE, ...)
      {
	  if(varcov && is.null(object$vcov))
	      object$vcov <- vcov(object, correlation=TRUE, sigm = object$sigma)
	  object
      })



### Plots for the ranef.mer class ----------------------------------------

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

if (FALSE) {
foo  <- function(x, data, ...)  ## old version of qqmath.ranef.mer
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
}## end{ unused }

qqmath.ranef.mer <- function(x, data, ...)
{
    prepanel.ci <- function(x, y, se, subscripts, ...) {
        x <- as.numeric(x)
        se <- as.numeric(se[subscripts])
        hw <- 1.96 * se
        list(xlim = range(x - hw, x + hw, finite = TRUE))
    }
    panel.ci <- function(x, y, se, subscripts, pch = 16, ...)  {
        panel.grid(h = -1,v = -1)
        panel.abline(v = 0)
        x <- as.numeric(x)
        y <- as.numeric(y)
        se <- as.numeric(se[subscripts])
        panel.segments(x - 1.96 * se, y, x + 1.96 * se, y, col = 'black')
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
            xyplot(rep(qnorm((rr - 0.5)/nr), ncol(x)) ~ unlist(x)[ord] | ind[ord],
                   se = se[ord], prepanel = prepanel.ci, panel = panel.ci,
                   scales = list(x = list(relation = "free")),
                   ylab = "Standard normal quantiles",
                   xlab = NULL, ...)
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
