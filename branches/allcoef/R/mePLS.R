##' Solve the penalized linear least squares problem associated with a
##' mixed-effects model.

##' @param Zt - transpose of sparse model matrix for the random effects
##' @param Lambda - relative covariance factor
##' @param y - response (or residual) vector
##' @param X - optional fixed-effects model matrix
##' @param L - optional sparse Cholesky factor
##' @param super - should a supernodal decomposition be used.  Ignored
##         if L is specified.

##' @return a list with a newly created or update CHMfactor object L
##' and the penalized least squares solution u.  If a fixed-effects
##' model matrix X is given then the factors RZX and RX are also
##' returned.

isLDL <- function(x)
{
    stopifnot(is(x, "CHMfactor"))
    as.logical(x@type[2])
}

mePLS <- function(Zt, Lambda, y, X, L, super = FALSE)
{
    stopifnot(is(Zt, "dgCMatrix"),
              is(Lambda, "sparseMatrix"))
    Ut <- crossprod(Lambda, Zt)
    if (missing(L)) {
        L <- Cholesky(tcrossprod(Ut), LDL = FALSE, super = super,
                      Imult = 1)
    } else {
        stopifnot(is(L, "CHMfactor"), isLDL(L))
        L <- update(L, Ut, 1)
        stopifnot(isLDL(L))
    }
    if (missing(X))
        return(list(L = L,
                    u = solve(L, Ut %*% as.numeric(y), system = "A")))
    PUtX <- solve(L, Ut %*% X, system = "P")
    RZX <- solve(L, PUtX, system = "L")
    RX <- chol(crossprod(X) - crossprod(RZX))
    list(L = L, RZX = RZX, RX = RX)
}
