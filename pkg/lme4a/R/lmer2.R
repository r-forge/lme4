##' Create a Matrix of 1's with prod(dim) = N and ncol = s
mkSqrtXwt <- function(N, s) {
    stopifnot((N %% s) == 0)
    Matrix(rep(1, N), nrow = N %/% s, ncol = s)
}

##' Create a reModule
##'
##' @param bars a list of parsed random-effects terms
##' @param fr a model frame in which to evaluate these terms
##'
##' @return a list of an reModule and a merTrms
mkReModule <- function(bars, fr, rwt = FALSE, s = 1L) {
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
    n <- nrow(fr)

    ## order terms stably by decreasing number of levels in the factor
    if (any(diff(nl)) > 0) {
        ord <- rev(order(nl))
        blist <- blist[ord]
        nl <- nl[ord]
    }
    Zt <- do.call(rBind, lapply(blist, "[[", "sm"))
    q <- nrow(Zt)
    ll <- list(Zt = Zt)

    ## Create and install Lambda, Lind, etc.  This must be done after
    ## any potential reordering of the terms.
    cnms <- lapply(blist, "[[", "cnms")
    nc <- sapply(cnms, length)          # no. of columns per term
    ncTrms <- nc
    nth <- as.integer((nc * (nc+1))/2)  # no. of parameters per term
    nb <- nc * nl                       # no. of random effects per term
    stopifnot(sum(nb) == q)

    ll$u <- numeric(q)
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
    ll$Class <- "reModule"
    if (rwt) {                          # different class with more information
        N <- ncol(Zt)
        ll$sqrtXwt <- mkSqrtXwt(N, s)
        if (s > 1) {
            ll$Ut <- Reduce(lapply(split(seq_len(N) , rep.int(1:s, rep.int(N %/% s, s))),
                                   function(cols) Ut[, cols]), "+")
        }
        ll$ubase <- ll$u
        ll$Class <- "reReMod"
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
    list(reMod = do.call("new", ll),
         Trms = new("merTrms", flist = fl, cnms = cnms))
} ## {mkReModule}

##' Create an feModule
##'
##' @param form Formula for the fixed-effects terms
##' @param fr a model frame in which to evaluate the formula
##' @param contrasts a list of constrasts for factors in fr
##' @param reMod the reModule for the model
##' @param sparseX Logical indicator of sparse X
##' @param rwt Logical indicator of reweightable
##' 
##' @return an object that inherits from feModule
mkFeModule <-
    function(form, fr, contrasts, reMod, sparseX, rwt = FALSE, s = 1L) {
                                        # fixed-effects model matrix X
    nb <- nobars(form[[3]])
    if (is.null(nb)) nb <- 1
    form[[3]] <- nb
    if (sparseX) {
        X <- sparse.model.matrix(form, fr, contrasts)
        ll <- list(Class = "lmerSpFeMod", X = X,
                   RX = Cholesky(crossprod(X), LDL = FALSE))
    } else {
        X <- Matrix(model.matrix(form, fr, contrasts))
        ll <- list(Class = "lmerDeFeMod", X = X, RX = chol(crossprod(X)))
    }
    rownames(ll$X) <- NULL
    ll$RZX <- reMod@Zt %*% ll$X
    ll$beta <- numeric(ncol(ll$X))
    if (rwt) {
        ll$V <- ll$X
        N <- nrow(ll$X)
        ll$sqrtXwt <- mkSqrtXwt(N, s)
        if (s > 1) {
            ll$V <- Reduce(lapply(split(seq_len(N) ,
                                        rep.int(1:s, rep.int(N %/% s, s))),
                                  function(rows) X[rows, ]), "+")
        }
        ll$betabase <- ll$beta
        ll$Class <- if (sparseX) "rwSpFeMod" else "rwDeFeMod"
    } else {                            # lmer model (the only non-reweightable type)
        ll$ZtX <- ll$RZX
        ll$XtX <- crossprod(ll$X)
    }
    do.call("new", ll)
}

mkRespMod <- function(fr, reMod, feMod, family = NULL, nlenv = NULL) {
    n <- nrow(fr)
                                        # components of the model frame
    y <- model.response(fr)
    # avoid problems with 1D arrays, but keep names
    if(length(dim(y)) == 1) {
        nm <- rownames(y)
        dim(y) <- NULL
        if(!is.null(nm)) names(y) <- nm
    }
    weights <- as.numeric(model.weights(fr))
    if (any(weights < 0))
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
    p <- ncol(feMod@X)
    q <- nrow(reMod@Zt)
    ll <- list(weights = weights, offset = offset, mu = numeric(n),
               resid = numeric(n), wrss = numeric(1), Utr = numeric(q),
               Vtr = numeric(p), cu = numeric(q), cbeta = numeric(p))
    if (is.null(family) && is.null(nlenv)) { # lmer
        ll$Class <- "lmerResp"
        y <- as.numeric(y)
        names(y) <- NULL
        ll$y <- y
        ll$Utr <- (reMod@Zt %*% y)@x
        ll$Vtr <- crossprod(feMod@X, y)@x
        return(do.call("new", ll))
    }
    if (!is.null(family)) {
        if (is.character(family))
            family <- get(family, mode = "function", envir = parent.frame(3))
        if (is.function(family)) family <- family()
        etastart <- unname(as.numeric(model.extract(fr, "etastart")))
        mustart <- unname(as.numeric(model.extract(fr, "mustart")))
        nobs <- n
        with(ll, eval(family$initialize))
    }
}


lmer2 <-
    function(formula, data, REML = TRUE, sparseX = FALSE,
             control = list(), start = NULL, verbose = 0, doFit = TRUE,
             compDev = TRUE, optimizer = c("nlminb", "bobyqa", "optimize"),
             subset, weights, na.action, offset,
             contrasts = NULL, ...)
{
    mf <- mc <- match.call()
### '...' handling up front, safe-guarding against typos ("familiy") :
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
    reList <- mkReModule(findbars(formula[[3]]), fr)
                                        # fixed-effects module
    feMod <- mkFeModule(formula, fr, contrasts, reList$reMod, sparseX)
    respMod <- mkRespMod(fr, reList$reMod, feMod)
    new("lmerTrms2", trms = reList$Trms, re = reList$reMod,
        fe = feMod, resp = respMod, REML = REML)
}

S4toEnv <- function(from) {
    stopifnot(isS4(from))
    ans <- new.env()
    for (nm in slotNames(from)) {
        sl <- slot(from, nm)
        for (nnm in slotNames(sl)) assign(nnm, slot(sl, nnm), envir = ans)
    }
    ans
}

setAs("lmer2", "optenv", function (from)
  {
      rho <- S4toEnv(from)
      rho$ldL2 <- numeric(1)
      rho$ldRX2 <- numeric(1)
      n <- length(rho$y)
      rho$nobs <- n
      rho$REML <- from@REML
      rho$nmp <- n - length(rho$beta)
      sP <- function(x) {
          stopifnot(length(x) == length(theta))
          ## new theta := x
          theta <<- as.numeric(x)
          Lambda@x[] <<- theta[Lind]           # update Lambda
          Ut <<- crossprod(Lambda, Zt)
          Matrix:::destructive_Chol_update(L, Ut, Imult = 1)
          cu <- solve(L, solve(L, crossprod(Lambda, Utr), sys = "P"), sys = "L")
          RZX <<- solve(L, solve(L, crossprod(Lambda, ZtX), sys = "P"), sys = "L")
          if (is(X, "sparseMatrix")) {
              RX <<- update(RX, XtX - crossprod(RZX))
              beta[] <<- solve(RX, Xty - crossprod(RZX, cu))
          } else {
              RX <<- chol(XtX - crossprod(RZX))
              beta[] <<- solve(RX, solve(t(RX), Vtr - crossprod(RZX, cu)))@x
          }
          u[] <<- solve(L, solve(L, cu - RZX %*% beta, sys = "Lt"), sys = "Pt")@x
          mu[] <<- (if(length(offset)) offset else 0) + (crossprod(Ut, u) + X %*% beta)@x
          pwrss <<- sum(c(y - mu, u)^2) # penalized residual sum of squares
          ldL2[] <<- 2 * determinant(L)$mod
          if (REML) {
              ldRX2[] <<- 2 * determinant(RX)$mod
              return(ldL2 + ldRX2 + nmp * (1 + log(2 * pi * pwrss/nmp )))
          }
          ldL2 + nobs * (1 + log(2 * pi * pwrss/nobs))
      }
      gP <- function() theta
      gB <- function() cbind(lower = lower,
                             upper = rep.int(Inf, length(theta)))
      environment(sP) <- environment(gP) <- environment(gB) <- rho
      new("optenv", setPars = sP, getPars = gP, getBounds = gB)
  })
