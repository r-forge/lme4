##' Determine if a CHMfactor object is LDL or LL
##' @param x - a CHMfactor object
##' @return TRUE if x is LDL, otherwise FALSE
isLDL <- function(x)
{
    stopifnot(is(x, "CHMfactor"))
    as.logical(x@type[2])
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
