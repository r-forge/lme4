## Class definitions for the package

setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

## TODO: export
setClass("lmList.confint", contains = "array")

## -------------------- lmer-related Classes --------------------------------

setOldClass("data.frame")
setOldClass("family")
setOldClass("logLik")


### Environment-based classes.  These will eventually replace the
### previous classes

##' Optimization environment class.

setClass("optenv", representation(setPars = "function",
                                  getPars = "function",
                                  getBounds = "function"))


##' Basic properties of a mixed-effects representation.
##'
##' The shared environment should contain objects y, X, Zt, Ut, fixef,
##' u, Lambda, Lind, theta, L and ldL2.
##'
setClass("merenv", representation("VIRTUAL"), contains = "optenv",
         validity = function(object)
     {
         rho <- env(object)
         if (!(is.numeric(y <- rho$y) &&
               (n <- length(y)) > 0))
             return("environment must contain a non-trivial numeric response y")
         if (!(is(X <- rho$X, "dMatrix") &&
               is(Zt <- rho$Zt, "dMatrix") &&
               ((N <- nrow(X)) == ncol(Zt))))
             return("environment must contain Matrix objects X and Zt with nrow(X) == ncol(Zt)")
         if (N %% n || N <= 0)
             return(sprintf("nrow(X) = %d must be a positive multiple of length(y) = %d",
                            N, n))
         p <- ncol(X)
         q <- nrow(Zt)
         if (!(is(Ut <- rho$Ut, "dMatrix") &&
               all(dim(Ut) == dim(Zt))))
             return("environment must contain Ut of same dimensions as Zt")
         if (!(is.numeric(fixef <- rho$fixef) && length(fixef) == p))
             return(sprintf("environment must contain a numeric vector fixef of length %d",
                            ncol(X)))
         if (!(is.numeric(u <- rho$u) &&
               length(u) == q))
             return(sprintf("environment must contain a numeric vector u of length %d",
                            q))
         if (!(is(Lambda <- rho$Lambda, "dMatrix") &&
               all(dim(Lambda) == q)))
             return(sprintf("environment must contain a %d by %d Matrix Lambda",
                            q, q))
         if (!(is.integer(Lind <- rho$Lind) &&
               length(Lind) == length(Lambda@x) &&
               min(Lind) == 1L))
             return(sprintf("environment must contain an integer vector Lind of length %d with minimum 1",
                            length(Lambda@x)))
         nth <- max(Lind)
         if (!(is.numeric(theta <- rho$theta) &&
               length(theta) == nth))
             return(sprintf("environment must contain a numeric vector theta of length %d",
                            nth))
         if (!(all(seq_along(theta) %in% Lind)))
             return("not all indices of theta occur in Lind")
         if (!(is(L <- rho$L, "CHMfactor") &&
               all(dim(L) == q)))
             return("environment must contain a CHMfactor L")
         TRUE
     })

##' Mixed-effects model representation based on random-effects terms
##'
##' The general merenv class does not associate components of the
##' random-effects vector, b, with particular terms in a formula.
##' In this class the random effects are associated with a set of
##' terms with grouping factors in the list flist.  The number of
##' columns in each term is available in the nc vector.
##'
##'
setClass("merenvtrms", representation("VIRTUAL"), contains = "merenv",
         validity = function(object) .Call(merenvtrms_validate, env(object)))
     ## {
     ##     rho <- env(object)
     ##     if (!(is.list(flist <- rho$flist) &&
     ##           all(sapply(flist, is.factor))))
     ##         return("environment must contain a list of factors, flist")
     ##     flseq <- seq_along(flist)
     ##     if (!(is.integer(asgn <- attr(flist, "assign")) &&
     ##           all(flseq %in% asgn) &&
     ##           all(asgn %in% flseq)))
     ##         return("asgn attribute of flist missing or malformed")
     ##     nl <- sapply(flist, function(x) length(levels(x)))[asgn]
     ##     if (!(is.list(cnms <- rho$cnms) &&
     ##           all(sapply(cnms, is.character)) &&
     ##           all(sapply(cnms, length) > 0) &&
     ##           length(cnms) == length(asgn)))
     ##         return("list of column names, cnms, must match asgn attribute in length")
     ## })

##' Linear mixed-effects model representation.
##'
##'
##'
##'
setClass("lmerenv", contains = "merenvtrms",
         validity = function(object)
     {
         rho <- env(object)
         p <- length(rho$fixef)
         q <- length(rho$u)
         if (!(is(Zty <- rho$Zty, "dMatrix") &&
               nrow(Zty) == q))
             return("environment must contain a column Matrix Zty")
         if (!(is(ZtX <- rho$ZtX, "dMatrix") &&
               all(dim(ZtX) == c(q, p))))
             return("environment must contain a q by p Matrix ZtX")
         if (!(is(RZX <- rho$RZX, "dMatrix") &&
               all(dim(RZX) == c(q, p))))
             return("environment must contain a q by p Matrix RZX")
         RX <- rho$RX
         if (is(RX, "CHMfactor")) RX <- as(RX, "sparseMatrix")
         if (!(is(RX, "dMatrix") && all(dim(RX) == c(p, p))))
             return("environment must contain a p by p Matrix RX")
         if (!(is.numeric(Xty <- rho$Xty) && length(Xty) == p))
             return("environment must contain a numeric p-vector Xty")
         if (!(is(XtX <- rho$XtX, "dMatrix") &&
               is(XtX, "symmetricMatrix") &&
               all(dim(RZX) == c(q, p))))
             return("environment must contain a symmetric Matrix XtX")
         TRUE
     })

setClass("glmerenv", contains = "merenvtrms",
         validity = function(object)
     {
         rho <- env(object)
         n <- length(rho$y)
         if (!(is.numeric(rho$mu) && length(rho$mu) == n))
             return("environment must contain a numeric vector \"mu\"")
         if (!(is.numeric(rho$muEta) && length(rho$muEta) == n))
             return("environment must contain a numeric vector \"muEta\"")
         if (!(is.numeric(rho$var) && length(rho$var) == n))
             return("environment must contain a numeric vector \"var\"")
         if (!(is.numeric(rho$sqrtrwt) && length(rho$sqrtrwt) == n))
             return("environment must contain a numeric vector \"sqrtrwt\"")
         if (!(is.list(family <- rho$family)))
             return("environment must contain a list \"family\"")
         if (!(is.character(family$family) && length(family$family) == 1))
             return("family list must contain a character variable \"family\"")
         if (!(is.character(family$link) && length(family$link) == 1))
             return("family list must contain a character variable \"link\"")
         TRUE
     })

