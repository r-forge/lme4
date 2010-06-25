## Class definitions for the package

setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

## TODO: export
setClass("lmList.confint", contains = "array")

## -------------------- lmer-related Classes --------------------------------

setOldClass("family")
## and  "data.frame', "logLik", "environment" are already defined

##' Random-effects module.
##'
##' Zt is the transpose of the sparse model matrix.  The number of
##' rows in Zt may be a multiple of the number of columns in Ut.

setClass("reModule",
         representation(L = "CHMfactor",
                        Lambda = "dgCMatrix",
                        Lind = "integer",
## I think that Ut can be removed except that it is still used in "simulate" and "bootMer"
                        Ut = "dgCMatrix",## U := Z Lambda; Ut := U' = Lambda' Z' = Lambda' Zt
                        Zt = "dgCMatrix",## = Z'
                        lower = "numeric",
                        theta = "numeric",
                        u = "numeric"),
         validity = function(object) {
             q <- nrow(object@Zt)
             if (!all(dim(object@Lambda) == q))
                 return("Lambda must be q by q where q = nrow(Zt)")
             if (nrow(object@Ut) != q || nrow(object@L) != q)
                 return("Number of rows in Zt, L and Ut must match")
             if (length(object@u) != q)
                 return("length(u) must be q = nrow(Zt)")
             if (length(object@Lind) != length(object@Lambda@x))
                 return("length(Lind) != length(Lambda@x)")
             if (!all(object@Lind %in% seq_along(object@theta)))
                 return("elements of Lind must be in 1:length(theta)")
             if (isLDL(object@L))
                 return("L must be an LL factor, not LDL")
             if (length(object@lower) != length(object@theta))
                 return("lengths of lower and theta must match")
             TRUE
         })

##' Random-effects module derived from random-effects terms
##'
##' In general an reModule does not associate components of the
##' random-effects vector, b, with particular terms in a formula. That
##' association is represented separately by this class.
##'
##'
setClass("reTrms",
         representation(flist = "list", cnms = "list"),
         contains = "reModule",
         validity = function(object)
     {
         flLen <- length(flist <- object@flist)
         if (flLen < 1 && !all(sapply(flist, is.factor)))
             return("flist must be a non-empty list of factors")
         l1 <- length(flist[[1]])
         if (!all(sapply(flist, function(el) length(el) == l1)))
             return("all factors in flist must have the same length")
         flseq <- seq_along(flist)
         if (!(is.integer(asgn <- attr(flist, "assign")) &&
               all(flseq %in% asgn) &&
               all(asgn %in% flseq)))
             return("asgn attribute of flist missing or malformed")
         if (!all(sapply(cnms <- object@cnms, is.character)) &&
             all(sapply(cnms, length) > 0) &&
             length(cnms) == length(asgn))
             return("list of column names, cnms, must match asgn attribute in length")
         nlev <- sapply(flist, function(fac) length(levels(fac)))
         nc <- sapply(object@cnms, length)
         q <- nrow(object@Zt)
         if (sum(nc * nlev[asgn]) != q)
             return("inconsistent dimensions in trms and re slots")
         TRUE
     })

##' Fixed-effects module
setClass("feModule", representation(beta = "numeric", "VIRTUAL"))

.feValid <- function(object) {
    p <- ncol(object@X)
    if (ncol(object@RZX) != p || ncol(object@RX) != p)
        return("Number of columns in X, RZX, and RX must match")
    TRUE
}

##' Dense fixed-effects module
setClass("deFeMod",
         representation(RZX = "dgeMatrix",
                        RX  =  "Cholesky",
			X   = "dgeMatrix"),
         contains = "feModule",
         validity = .feValid)

##' Sparse fixed-effects module
setClass("spFeMod",
         representation(RZX = "dgCMatrix",
                        RX =  "CHMfactor",
                        X =   "dgCMatrix"),
         contains = "feModule",
         validity = function(object) {
             if (is(rr <- .feValid(object), "character")) return(rr)
             if (isLDL(object@RX))
                 return("RX must be an LL factor, not LDL")
             TRUE
         })

##' mer response module
##' y, offset and mu are as expected.  Note that length(offset) can be a multiple of length(y)
##' weights are the prior weights
##' sqrtrwt and sqrtXwt are the square roots of residual and X weights
##' wtres is the vector of weighted residuals
setClass("merResp",
         representation(mu = "numeric",
                        offset = "numeric",
                        sqrtXwt = "matrix",
                        sqrtrwt = "numeric", # sqrt(residual weights)
                        weights = "numeric", # prior weights
                        y = "numeric"),
         validity = function(object) {
             n <- length(object@y)
             if (any(n != sapply(lapply(c("weights","sqrtrwt","mu"#,"wtres"
                     ), slot, object = object), length)))
                 return("lengths of weights, sqrtwt and mu must match length(y)")
             lo <- length(object@offset)
             if (!lo || lo %% n)
                 return("length(offset) must be a positive multiple of length(y)")
             if (length(object@sqrtXwt) != lo)
                 return("length(sqrtXwt) must equal length(offset)")
             if (nrow(object@sqrtXwt) != n)
                 return("nrow(sqrtXwt) != length(y)")
             TRUE
         })

setClass("lmerResp", representation(REML = "integer"),
         contains = "merResp")

##' glmer response module
setClass("glmerResp",
         representation(family =  "family",
                        eta =    "numeric",
                        n =      "numeric"), # for evaluation of the aic
         contains = "merResp",
         validity = function(object) {
             if (length(object@eta) != length(object@y))
                 return("lengths of eta and y must match")
         })

##' nlmer response module
setClass("nlmerResp",
         representation(nlenv = "environment",
                        nlmod = "call",
                        pnames = "character"),
         contains = "merResp",
         validity = function(object) {
             n <- length(object@y)
             N <- length(object@offset)
             s <- N %/% n
             lpn <- length(object@pnames)
             if (lpn != s) return(sprintf("length(pnames) = %d != s = %d", lpn, s))
             dd <- dim(object@sqrtXwt)
             if (!all(dd == c(n, s))) {
                 return(sprintf("dim(gradient) = (%d, %d), n = %d, s = %d",
                                dd[1], dd[2], n, s))
             }
             TRUE
         })

##' nglmer response module
setClass("nglmerResp", contains = c("glmerResp", "nlmerResp"))

setClass("merMod",
         representation(call    = "call",
                        devcomp = "list",
			frame   = "data.frame", # "model.frame" is not S4-ized yet
                        re      = "reModule",
                        fe      = "feModule",
                        resp    = "merResp"))
