##' Determine if a CHMfactor object is LDL or LL
##' @param x - a CHMfactor object
##' @return TRUE if x is LDL, otherwise FALSE
isLDL <- function(x)
{
    stopifnot(is(x, "CHMfactor"))
    as.logical(x@type[2])
}


lmer2 <-
    function(formula, data, family = gaussian, REML = TRUE,
             control = list(), start = NULL, verbose = FALSE,
             subset, weights, na.action, offset, contrasts = NULL,
             model = TRUE, mustart, etastart, ...)
{
    if (missing(data)) data <- environment(formula)
    stopifnot(length(formula <- as.formula(formula)) == 3)
                                        # evaluate and install the model frame
    mf <- mc <- match.call()
    m <- match(c("data", "subset", "weights", "na.action",
                 "etastart", "mustart", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula)
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    n <- nrow(fr)
                                        # components of the model frame
    y <- model.response(fr)
    if(length(dim(y)) == 1) { # avoid problems with 1D arrays, but keep names
        nm <- rownames(y)
        dim(y) <- NULL
        if(!is.null(nm)) names(y) <- nm
    }
    weights <- as.numeric(as.vector(model.weights(fr)))
    if (length(weights) == 0)
        weights <- rep.int(1, n)
    if (any(weights < 0))
        stop(gettext("negative weights not allowed", domain = "R-lme4"))
    loff <- length(offset <- as.numeric(as.vector(model.offset(fr))))
    if (loff) {
        if (loff == 1) {
            offset <- rep.int(offset, n)
        } else if (loff != nobs) {
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                          loff, n), domain = "R-lme4")
        }
    }
    ##            # starting values expressed as mu or eta (glmer only)
    ## mustart <- model.extract(mf, "mustart")
    ## etastart <- model.extract(mf, "etastart")
    
    fe.form <- formula           # evaluate fixed-effects model matrix
    nb <- nobars(formula[[3]])   # fixed-effects terms only
    if (is.null(nb)) nb <- 1
    fe.form[[3]] <- nb
    X <- Matrix(model.matrix(fe.form, fr, contrasts))
    rownames(X) <- NULL
    
    p <- ncol(X)
    stopifnot((nmp <- n - p) > 0)
    beta <- numeric(p)
    names(beta) <- colnames(X)
                                        # enforce modes on some vectors
    y <- unname(as.double(y))   # must be done after initialize
    attr(frame, "terms") <- NULL
    ## mustart <- unname(as.double(mustart))
    ## etastart <- unname(as.double(etastart))    
    ## if (exists("n", envir = rho))
    ##     n <- as.double(n)
    ## if (family$family %in% c("binomial", "poisson"))
    ##     dims["useSc"] <- 0L
    
    ## Check for method argument which is no longer used
    if (!is.null(method <- list(...)$method)) {
        msg <- paste("Argument", sQuote("method"),
                     "is deprecated.\nUse", sQuote("nAGQ"),
                     "to choose AGQ.  PQL is not available.")
        if (match.arg(method, c("Laplace", "AGQ")) == "Laplace") {
            warning(msg)
        } else stop(msg)
    }
    
    bars <- expandSlash(findbars(formula[[3]]))
    if (!length(bars)) stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))

    stopifnot(all(sapply(bars, "[[", 2) == 1)) # check for simple, scalar terms
    flist <- lapply(bars,
                    function(x)
                    eval(substitute(factor(fac), list(fac = x[[3]])), fr))
    rm(nb, mf, fe.form, loff, bars)
    
    RX <- chol(XtX <- crossprod(X))     # check for full column rank
    Xty <- unname(as.vector(crossprod(X, y)))
                  
    Ut <- Zt <- do.call(rBind, lapply(flist, as, "sparseMatrix"))
    RZX <- Ut %*% X
    Sind <- rep.int(seq_along(flist),
                    sapply(flist, function(x) length(levels(x))))
    u <- numeric(nrow(Zt))
    theta <- numeric(length(flist))
    fitted <- y
    prss <- 0
    ldL2 <- 0
    
    L <- Cholesky(tcrossprod(Zt), LDL = FALSE, Imult = 1)
    S <- Diagonal(x = theta[Sind])
    new("merenv",
        setPars = function(x)
     {
         theta <<- as.numeric(x)
         stopifnot(length(theta) == length(flist))
         S@x[] <- theta[Sind]           # update S
         Ut <<- crossprod(S, Zt)
         L <<- update(L, Ut, mult = 1)
         cu <- solve(L, solve(L, Ut %*% y, sys = "P"), sys = "L")
         RZX <<- solve(L, solve(L, Ut %*% X, sys = "P"), sys = "L")
         RX <<- chol(XtX - crossprod(RZX))
         cb <- solve(t(RX), Xty - crossprod(RZX, cu))
         beta[] <<- as.vector(solve(RX, cb))
         u[] <<- as.vector(solve(L, solve(L, cu - RZX %*% beta, sys = "Lt"), sys = "Pt"))
         fitted[] <<- as.vector(crossprod(Ut, u) + X %*% beta)
         prss <<- sum(c(y - fitted, as.vector(u))^2) # penalized residual sum of squares
         ldL2 <<- as.vector(determinant(L)$mod)
         if (REML) return(as.vector(ldL2 + 2*determinant(RX)$mod +
                                    nmp * (1 + log(2 * pi * prss/nmp))))
         ldL2 + n * (1 + log(2 * pi * prss/n))
     },
         getPars = function() theta,
         getBounds = function() list(lower = rep.int(0, length(theta)), upper = rep.int(Inf, length(theta)))
         )
}

##' Return a function to evaluate the profiled deviance or REML
##' criterion for a linear mixed model with simple, scalar random
##' effects terms

##' @param flist a list of factors from which to generate simple,
##'              scalar random effects terms
##' @param y the numeric response vector
##' @param X model matrix for the fixed-effects parameters
##' @param REML optional logical scalar indicating if the REML
##'             criterion should be used.  (default is TRUE).
##' @param super optional logical scalar indicating if a supernodal
##'        decomposition should be used (default is FALSE).

##' @return a function closure that evaluates the criterion chosen.
##'   Calling optimize or nlminb with this function as an argument
##'   optimizes the criterion.  The updated parameters and derived
##'   quantities are in the enclosing environment for this function.
simplemer <- function(flist, y, X, REML = TRUE, super = FALSE)
{
    X <- Matrix(X)
    n <- length(y <- as.numeric(y))

    stopifnot(n > 0,                 # check arguments for consistency
              nrow(X) == n,
              is.list(flist),
              length(flist) > 0,
              all(sapply(flist, is.factor)),
              all(sapply(flist, length) == n))
    super <- as.logical(super)[1]
    REML <- as.logical(REML)[1]
    nmp <- n - (p <- ncol(X))
    beta <- numeric(p)

    RX <- chol(XtX <- crossprod(X))     # check for full column rank
    Xty <- crossprod(X, y)
    
    Ut <- Zt <- do.call(rBind, lapply(flist, as, "sparseMatrix"))
    RZX <- Ut %*% X
    Sind <- rep.int(seq_along(flist),
                     sapply(flist, function(x) length(levels(factor(x)))))
    u <- numeric(nrow(Zt))
    theta <- numeric(length(flist))
    fitted <- y
    prss <- 0
    ldL2 <- 0
    
    L <- Cholesky(tcrossprod(Zt), LDL = FALSE, Imult = 1, super = super)
    S <- Diagonal(x = theta[Sind])
    new("merenv",
        setPars = function(x)
     {
         theta <<- as.numeric(x)
         stopifnot(length(theta) == length(flist))
         S@x[] <- theta[Sind]           # update S
         Ut <<- crossprod(S, Zt)
         L <<- update(L, Ut, mult = 1)
         cu <- solve(L, solve(L, Ut %*% y, sys = "P"), sys = "L")
         RZX <<- solve(L, solve(L, Ut %*% X, sys = "P"), sys = "L")
         RX <<- chol(XtX - crossprod(RZX))
         cb <- solve(t(RX), Xty - crossprod(RZX, cu))
         beta <<- solve(RX, cb)
         u <<- solve(L, solve(L, cu - RZX %*% beta, sys = "Lt"), sys = "Pt")
         fitted <<- as.vector(crossprod(Ut, u) + X %*% beta)
         prss <<- sum(c(y - fitted, as.vector(u))^2) # penalized residual sum of squares
         ldL2 <<- as.vector(determinant(L)$mod)
         if (REML) return(as.vector(ldL2 + 2*determinant(RX)$mod +
                                    nmp * (1 + log(2 * pi * prss/nmp))))
         ldL2 + n * (1 + log(2 * pi * prss/n))
     },
         getPars = function() theta,
         getBounds = function() list(lower = rep.int(0, length(theta)), upper = rep.int(Inf, length(theta)))
         )
}

##' Solve the penalized linear least squares problem associated with a
##' mixed-effects model.

##' @param Ut - crossprod(Lambda, Zt) where Zt is the transpose of
##' sparse model matrix for the random effects and Lambda is the
##' left relative covariance factor
##' @param r - response (or residual) vector
##' @param V - optional fixed-effects model matrix
##' @param L - optional sparse Cholesky factor compatible with Ut
##' @param u0 - optional current value of u
##' @param super - logical scalar, should a supernodal decomposition
##' be used? Ignored if L is specified.

##' @return a list with a newly created or update CHMfactor object L
##' and the penalized least squares solution u.  If a fixed-effects
##' model matrix X is given then the factors RZX and RX are also
##' returned.

mePLS <- function(Ut, r, V, L, u0, super = FALSE)
{
    stopifnot(is(Ut, "dgCMatrix"))
    if (missing(L)) {
        L <- Cholesky(tcrossprod(Ut), LDL = FALSE, super = super,
                      Imult = 1)
    } else {
        stopifnot(is(L, "CHMfactor"), isLDL(L))
        L <- update(L, Ut, 1)
        stopifnot(isLDL(L))
    }
    Utr <- Ut %*% as.numeric(r)
    if (!missing(u0)) Utr <- Utr - u0
    if (missing(V))
        return(list(L = L,
                    u = solve(L, Utr, system = "A")))
    cu <- solve(L, solve(L, Utr, system = "P"), system = "L")
    RZX <- solve(L, solve(L, Ut %*% V, system = "P"), system = "L")
    ## Martin: is there a better way to do this for the general case?
    ## That is, should we allow a downdate operation for a Cholesky
    ## factor?
    RX <- chol(crossprod(V) - crossprod(RZX))
    cb <- solve(t(RX), crossprod(V, r) - crossprod(RZX, cu))
    delb <- solve(RX, cb)
    delu <- solve(L, solve(L, cu - RZX %*% delb, system = "Pt"), system = "Lt")
    list(L = L, RZX = RZX, RX = RX, cu = cu, cb = cb, delu = delu,
         delb = delb)
}

## Generate a model with one or more simple scalar random effects
## terms
genSimple <- function(flist, y, X)
{
    if (is.factor(flist)) flist <- list(flist)
    n <- length(y <- as.numeric(y))
    flist <- lapply(flist, factor)      # drops unused levels
    stopifnot(all(sapply(flist, length) == n),
              is.matrix(X) || is(X, "Matrix"),
              nrow(X) == n)
    k <- length(flist)
    theta <- rep.int(1, k)
    nlev <- sapply(flist, function(f) length(levels(f)))
    q <- sum(nlev)
    Zt <- drop0(do.call(rBind, lapply(flist, as, Class =
                                      "sparseMatrix")))
    Lambda <- Diagonal(q)
    Ut <- crossprod(Lambda, Zt)
    L <- mePLS(Ut, y)$L
    getPars <- function() theta
    getBounds <- function()
        list(lower = rep.int(0, k), upper = rep.int(Inf, k))
    setPars <- function(th)
    {
        theta <<- th
        stopifnot(is.numeric(th), length(th) == k)
        Lambda <<- Diagonal(x = rep.int(th, nlev))
        Ut <<- crossprod(Lambda, Zt)
        ll <- mePLS(Ut, y, X, L)
        fitd <- crossprod(Ut, ll$delu) + X %*% ll$delb
        rsqr <- sum(c(y - as.vector(fitd), as.vector(ll$delu))^2)
        as.vector(determinant(ll$L, logarithm = TRUE)$modulus) +
            n * (1 + log(2*pi*rsqr/n))
    }
    list(getPars = getPars, getBounds = getBounds, setPars = setPars)
}        
        
simpleDev <- function(Zt, y, X)
{
    y <- as.numeric(y)
    L <- mePLS(Zt, y)$L                 # initialize the decomposition
    n <- ncol(Zt)
    stopifnot(is.matrix(X) || is(X, "Matrix"),
              nrow(X) == n,
              length(y) == n)
    devfunc <- function(theta)
    {
        Ut <- as.numeric(theta)[1] * Zt
        ll <- mePLS(Ut, y, X, L)
        fitd <- crossprod(Ut, ll$delu) + X %*% ll$delb
        rsqr <- sum(c(y - as.vector(fitd), as.vector(ll$delu))^2)
        determinant(ll$L, logarithm = TRUE)$modulus +
            n * (1 + log(2*pi*rsqr/n))
    }
    optimize(devfunc, c(0, 5))
}

## An example would be
# data(Dyestuff, package = "lme4")
# Zt <- as(Dyestuff$Batch, "sparseMatrix")
# y <- Dyestuff$Yield
# X <- matrix(rep.int(1, length(y)), nc = 1)
# simpleDev(Zt, y, X)
#
## except that doesn't seem to give the expected answers.  The
# deviance calculations match those in the book at zero and at one but
# not in between 0 and 1.
