.f <- if(package_version(packageDescription("Matrix")$Version) >=
         "0.999375-31") 2 else 1
## Also suppress the warning
assign("det_CHMfactor.warn", TRUE, envir = Matrix:::.MatrixEnv)

##' Install various objects, especially Lambda, Lind, Zt and Ut in
##' environment rho based on the list bars of random effects terms and
##' model frame fr.
##' @param bars a list of parsed random-effects terms
##' @param fr a model frame in which to evaluate these terms
##' @rho an environment that is modified by this function
makeZt <- function(bars, fr, rho)
{
    stopifnot(is.list(bars), all(sapply(bars, is.language)),
              inherits(fr, "data.frame"), is.environment(rho))
    if (!length(bars))
        stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    blist <- lapply(bars,
                    function(x)
                {
                    ff <- eval(substitute(factor(fac), list(fac = x[[3]])), fr)
                    nl <- length(levels(ff))
                    mm <- model.matrix(eval(substitute( ~ foo,
                                                       list(foo = x[[2]]))), fr)
                    nc <- ncol(mm)
                    nseq <- seq_len(nc)
                    sm <- as(ff, "sparseMatrix")
                    if (nc  > 1) 
                        sm <- do.call(rBind, lapply(nseq, function(i) sm))
                    sm@x[] <- t(mm[])
                    list(ff = ff, sm = sm, nc = nc, nl = nl)
                })
    ## number of random effects for each term
    nlev <- sapply(blist, function(el) length(levels(el$ff)))
    ## order terms stably by decreasing number of levels in the factor
    if (any(diff(nlev)) > 0) blist <- blist[rev(order(nlev))]
                                        # massage the factor list
    fl <- lapply(blist, "[[", "ff")
    asgn <- seq_along(fl)
                                        # check for repeated factors
    fnms <- names(fl)
    if (length(fnms) > length(ufn <- unique(fnms))) {
### FIXME: check that the lengths of the number of levels coincide
        fl <- fl[match(ufn, fnms)]
        asgn <- match(fnms, ufn)
    }
    names(fl) <- ufn
    fl <- do.call(data.frame, c(fl, check.names = FALSE))
    attr(fl, "assign") <- asgn
    
    rho$flist <- fl
    rho$Zt <- do.call(rBind, lapply(blist, "[[", "sm"))
    rho$Ut <- rho$Zt
    nt <- length(blist)                 # no. of r.e. terms
    nl <- sapply(blist, "[[", "nl")     # no. of levels per term
    nc <- sapply(blist, "[[", "nc")     # no. of columns per term
    nth <- as.integer((nc * (nc+1))/2)  # no. of parameters per term
    rho$q <- q <- sum(nb <- nc * nl)    # total no. of random effects
    rho$u <- numeric(q)
    rho$theta <- numeric(sum(nth))
    boff <- cumsum(c(0L, nb))           # offsets into b
    thoff <- cumsum(c(0L, nth))         # offsets into theta
    lst <- lapply(seq_along(blist), function(i)
              {
                  n <- nc[i] * nl[i]
                  mm <- matrix(seq_len(n), nc = nc[i])
                  dd <- diag(nc[i])
                  ltri <- lower.tri(dd, diag = TRUE)
                  ii <- row(dd)[ltri]
                  jj <- col(dd)[ltri]
                  dd[cbind(ii, jj)] <- seq_along(ii)
                  list(i = as.vector(mm[, ii]) + boff[i],
                       j = as.vector(mm[, jj]) + boff[i],
                       x = rep.int(seq_along(ii), rep.int(nl[i], length(ii))) +
                       thoff[i])
              })
    x <- unlist(lapply(lst, "[[", "x"))
    if (all(nc == 1L)) {
        rho$Lambda <- Diagonal(x = x)
    } else {
        rho$Lambda <-sparseMatrix(i = unlist(lapply(lst, "[[", "i")),
                                  j = unlist(lapply(lst, "[[", "j")),
                                  x = x)
    }
    rho$Lind <- as.integer(rho$Lambda@x)
    lower <- -Inf * (rho$theta + 1)
    lower[unique(diag(rho$Lambda))] <- 0
    rho$lower <- lower
}

lmer2 <-
    function(formula, data, family = gaussian, REML = TRUE, sparseX = FALSE,
             control = list(), start = NULL, verbose = FALSE,
             subset, weights, na.action, offset, contrasts = NULL, ...)
{
    mf <- mc <- match.call()
    if (!missing(family)) {      # call glmer if family is not missing
        mc[[1]] <- "glmer"
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
    rho$frame <- fr <- eval(mf, parent.frame())
    rho$n <- n <- nrow(fr)
                                        # components of the model frame
    y <- model.response(fr)
    if(length(dim(y)) == 1) { # avoid problems with 1D arrays, but keep names
        nm <- rownames(y)
        dim(y) <- NULL
        if(!is.null(nm)) names(y) <- nm
    }
    rho$y <- y <- unname(as.double(y))

    weights <- as.numeric(as.vector(model.weights(fr)))
    if (any(weights < 0))
        stop(gettext("negative weights not allowed", domain = "R-lme4"))
    rho$weights <- weights

    loff <- length(offset <- as.numeric((model.offset(fr))))
    if (loff) {
        if (loff == 1) {
            offset <- rep.int(offset, n)
        } else if (loff != n) {
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                          loff, n), domain = "R-lme4")
        }
    }
    rho$offset <- offset
    
    fe.form <- formula           # evaluate fixed-effects model matrix
    nb <- nobars(formula[[3]])   # fixed-effects terms only
    if (is.null(nb)) nb <- 1
    fe.form[[3]] <- nb
    X <- if (sparseX) {
        sparse.model.matrix(fe.form, fr, contrasts)
    } else Matrix(model.matrix(fe.form, fr, contrasts))
    rownames(X) <- NULL
    rho$X <- X
    
    rho$p <- p <- ncol(X)
    stopifnot((rho$nmp <- n - p) > 0)
    beta <- numeric(p)
    names(beta) <- colnames(X)
    rho$beta <- beta
    
    attr(fr, "terms") <- NULL
    
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
    
    rho$RX <- chol(rho$XtX <- crossprod(X)) # check for full column rank
    rho$Xty <- unname(as.vector(crossprod(X, y)))
    
    rho$RZX <- rho$Ut %*% X
    rho$fitted <- y
    rho$prss <- 0
    rho$ldL2 <- 0
    rho$ldRX2 <- 0
    
    rho$L <- Cholesky(tcrossprod(rho$Zt), LDL = FALSE, Imult = 1)
    rho$Zty <- rho$Zt %*% y
    rho$ZtX <- rho$Zt %*% X
    sP <- function(x)
    {
### FIXME: weights are not yet incorporated
        stopifnot(length(x) == length(theta))
        theta <<- as.numeric(x)
        Lambda@x[] <<- theta[Lind]           # update S
        Ut <<- crossprod(Lambda, Zt)
        Matrix:::destructive_Chol_update(L, Ut, Imult = 1)
        cu <- solve(L, solve(L, crossprod(Lambda, Zty), sys = "P"),
                    sys = "L")
        RZX <<- solve(L, solve(L, crossprod(Lambda, ZtX), sys = "P"),
                      sys = "L")
        RX <<- chol(XtX - crossprod(RZX))
        cb <- solve(t(RX), Xty - crossprod(RZX, cu))
        beta[] <<- solve(RX, cb)@x
        u[] <<- solve(L, solve(L, cu - RZX %*% beta, sys = "Lt"),
                      sys = "Pt")@x
        fitted[] <<- (if(length(offset)) offset else 0) +
            (crossprod(Ut, u) + X %*% beta)@x
        prss <<- sum(c(y - fitted, u)^2) # penalized residual sum of squares
        ldL2[] <<- .f * determinant(L)$mod
        ldRX2[] <<- 2 * determinant(RX)$mod
        if (!REML) return(ldL2 + n * (1 + log(2 * pi * prss/n)))
        ldL2 + ldRX2 + nmp * (1 + log(2 * pi * prss/nmp))
    }
    gP <- function() theta
    gB <- function() cbind(lower = rep.int(0, length(theta)),
                           upper = rep.int(Inf, length(theta)))
    environment(sP) <- environment(gP) <- environment(gB) <- rho
    new("merenv", setPars = sP, getPars = gP, getBounds = gB)
}

setMethod("fixef", "merenv", function(object, ...) env(object)$beta)

setMethod("ranef", "merenv", function(object, ...)
          {
              with(env(object), {
                  ans <- S@x * u
                  if (length(Tind)) ans <- (T %*% ans)@x
                  ans <- split(ans, Sind)
                  names(ans) <- names(flist)
                  ans
              })
          })

devcomp <- function(x, theta, ...)
{
    stopifnot(is(x, "merenv"))
    if (!missing(theta)) x@setPars(theta)
    with(env(x),
         list(cmp = c(ldL2 = ldL2, ldRX2 = ldRX2, prss = prss,
              deviance = ldL2 + n * (1 + log(2 * pi * prss/n)),
              REML = ldL2 + ldRX2 + nmp * (1 + log(2 * pi * prss/nmp))),
              dims = c(n = n, p = p, nmp = nmp, q = nrow(Zt))))
}

