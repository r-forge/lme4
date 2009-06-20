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

##' @param Zt - crossprod(Lambda, Zt) where Zt is the transpose of
##' sparse model matrix for the random effects and Lambda is the
##' left relative covariance factor
##' @param y - response (or residual) vector
##' @param V - optional fixed-effects model matrix
##' @param L - optional sparse Cholesky factor compatible with 
##' @param super - logical scalar, should a supernodal decomposition
##' be used? Ignored if L is specified.

##' @return a list with a newly created or update CHMfactor object L
##' and the penalized least squares solution u.  If a fixed-effects
##' model matrix X is given then the factors RZX and RX are also
##' returned.

mePLS <- function(Ut, y, V, L, super = FALSE)
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
    if (missing(V))
        return(list(L = L,
                    u = solve(L, Ut %*% as.numeric(y), system = "A")))
    PUtX <- solve(L, Ut %*% X, system = "P")
    RZX <- solve(L, PUtX, system = "L")
    RX <- chol(crossprod(X) - crossprod(RZX))
    list(L = L, RZX = RZX, RX = RX)
}
