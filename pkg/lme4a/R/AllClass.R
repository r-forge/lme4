## Class definitions for the package

setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

## TODO: export
setClass("lmList.confint", contains = "array")

## -------------------- lmer-related Classes --------------------------------

## setOldClass("family")
## and  "data.frame', "logLik", "environment" are already defined

##' Random-effects module.
##'
##' Zt is the transpose of the sparse model matrix.  The number of
##' columns in Zt may be a multiple of the number of columns in Ut.

setClass("reModule",
         representation(L = "CHMfactor",
                        Lambda = "dgCMatrix",
                        Lind = "integer",
## I think that Ut can be removed except that it is still used in "simulate" and "bootMer"
##                        Ut = "dgCMatrix",## U := Z Lambda; Ut := U' = Lambda' Z' = Lambda' Zt
                        Zt = "dgCMatrix",## = Z'
                        lower = "numeric",
                        theta = "numeric",
                        u = "numeric"),
         validity = function(object) {
             q <- nrow(object@Zt)
             if (!all(dim(object@Lambda) == q))
                 return("Lambda must be q by q where q = nrow(Zt)")
#             if (nrow(object@Ut) != q || nrow(object@L) != q)
#                 return("Number of rows in Zt, L and Ut must match")
             if (nrow(object@L) != q)                 
                 return("Number of rows in Zt and L must match")
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
         representation(flist = "list", cnms = "list"
##-nL                   , nLevs = "integer"
                        ),
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
##-nL 	 if(!all(nlev == object@nLevs))
##-nL 	     return("'nLevs' is not consistent length(levels(.)) of 'flist'")
         nc <- sapply(object@cnms, length)
         q <- nrow(object@Zt)
         if (sum(nc * nlev[asgn]) != q)
             return("inconsistent dimensions in trms and re slots")
         TRUE
     })

##' Fixed-effects module
#setClass("feModule", representation(coef = "numeric", "VIRTUAL"))

.feValid <- function(object) {
    p <- ncol(object@X)
    if (ncol(object@RZX) != p || ncol(object@fac) != p)
        return("Number of columns in X, RZX, and fac must match")
    TRUE
}

##' Dense fixed-effects module
setClass("deFeMod",
         representation(RZX = "dgeMatrix"),
##                      fac = "Cholesky",
##			X   = "dgeMatrix"),
         contains = "dPredModule",
         validity = .feValid)

##' Sparse fixed-effects module
setClass("spFeMod",
         representation(RZX = "dgCMatrix"),
##                      fac = "CHMfactor",
##                      X   = "dgCMatrix"),
         contains = "sPredModule",
         validity = function(object) {
             if (is(rr <- .feValid(object), "character")) return(rr)
             if (isLDL(object@fac))
                 return("fac must be an LL factor, not LDL")
             TRUE
         })

setClass("merMod",
         representation(call    = "call",
                        devcomp = "list",
			frame   = "data.frame", # "model.frame" is not S4-ized yet
                        re      = "reModule",
                        fe      = "predModule",
                        resp    = "respModule"))

setClass("lmerResp", representation(REML = "integer"), contains = "respModule")

setAs("Rcpp_reModule", "reModule", function(from)
  {
      new("reModule", L=from$L, Lambda=from$Lambda, Lind=from$Lind,
          Zt=from$Zt, lower=from$lower, theta=from$theta, u=from$u)
  })

setAs("Rcpp_deFeMod", "deFeMod", function(from)
  {
      new("deFeMod", RZX=from$RZX, X=from$X, fac=from$RX,
          coef=from$coef, Vtr=from$Vtr)
  })

setAs("Rcpp_lmerResp", "lmerResp", function(from)
  {
      new("lmerResp", REML=from$REML, mu=from$mu, offset=from$offset,
          sqrtXwt=from$sqrtXwt, sqrtrwt=from$sqrtrwt, weights=from$weights,
          wtres=from$wtres, y=from$y)
  })

