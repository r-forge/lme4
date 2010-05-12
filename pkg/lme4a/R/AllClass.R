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

### Environment-based classes.  These will eventually replace the
### previous classes

##' Optimization environment class.

setClass("optenv", representation(setPars = "function",
				  getPars = "function",
				  getBounds = "function"),
	 contains = "environment")


##' Basic properties of a mixed-effects representation.
##'
##' The shared environment should contain objects y, X, Zt, Ut, beta,
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
         if (!(is.numeric(beta <- rho$beta) && length(beta) == p))
             return(sprintf("environment must contain a numeric vector beta of length %d",
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
         p <- length(rho$beta)
         q <- length(rho$u)
	 if (!(is.numeric(Utr <- rho$Utr) && length(Utr) == q))
	     return("environment must contain a numeric q-vector Utr")
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
         if (!(is.numeric(Vtr <- rho$Vtr) && length(Vtr) == p))
             return("environment must contain a numeric p-vector Vtr")
         if (!(is(XtX <- rho$XtX, "dMatrix") &&
               is(XtX, "symmetricMatrix") &&
               all(dim(RZX) == c(q, p))))
             return("environment must contain a symmetric Matrix XtX")
         TRUE
     })

## Either a dense _or_ sparse Cholesky factorization
setClassUnion("CholKind",
              members = c("Cholesky", "CHMfactor"))

if(is.null(getClassDef("dsymmetricMatrix"))) ## not very latest version of 'Matrix':
    ## Virtual class of numeric symmetric matrices -- new (2010-04) - for lme4
    setClass("dsymmetricMatrix", representation("VIRTUAL"),
             contains = c("dMatrix", "symmetricMatrix"))

## Some __NON-exported__ classes just for auto-validity checking:
setClassUnion("dgC_or_diMatrix", members = c("dgCMatrix", "ddiMatrix"))
setClassUnion("dsC_or_dpoMatrix", members = c("dsCMatrix", "dpoMatrix"))

### FIXME?  This is currently modelled after "old-lme4"
### ----- following  g?l?merenv, should rather have *common* slots here,
##  and also define "lmer" , "glmer" etc  *classes* with the specific slots?
## __ FIXME __
setClass("mer",
	 representation(
			env = "merenv",
			frame = "data.frame",# model frame (or empty frame)
			call = "call",	 # matched call
			sparseX = "logical",
			X = "dMatrix", # fixed effects model matrix (sparse or dense)
			## Xst = "dgCMatrix", # sparse fixed effects model matrix
			Lambda = "dgC_or_diMatrix", # ddi- or dgC-Matrix"
			Zt = "dgCMatrix",# sparse form of Z'
			ZtX = "dMatrix", # [dge- or dgC-Matrix]
			Utr = "numeric", # Zty was "dgeMatrix",
			pWt = "numeric",# prior weights,   __ FIXME? __
			offset = "numeric", # length 0 -> no offset
			y = "numeric",	 # response vector
                                        # Gp = "integer",  # pointers to row groups of Zt

			devcomp = "list", ## list(cmp = ...,  dims = ...)
                        ## this *replaces* (?) previous 'dims' and 'deviance'
			##-- deviance = "numeric", # ML and REML deviance and components
			dims = "integer",# dimensions and indicators
			##-- dims : old-lme4's lmer(): named vector of length 18, names :
                        ## --> ../../lme4/R/lmer.R 'dimsNames' & 'dimsDefaults'
                        ## c("nt", "n", "p", "q", "s", "np", "LMM", "REML",
                        ##   "fTyp", "lTyp", "vTyp", "nest", "useSc", "nAGQ",
                        ##   "verb", "mxit", "mxfn", "cvg")

			##--- slots that vary during optimization------------

			V = "matrix",	 # gradient matrix
			Ut = "dgCMatrix", #
			L = "CHMfactor", # Cholesky factor of weighted P(.. + I)P'
			Lind = "integer",
			beta = "numeric" ,# fixed effects (length p)
			ranef = "numeric",# random effects (length q)
			u = "numeric",	 # orthogonal random effects (q)
			theta = "numeric",# parameter for Sigma or Lambda
			mu = "numeric",	 # fitted values at current beta and b
			gamma = "numeric",#  for NLME - otherwise == mu (?)
			##? resid = "numeric",# raw residuals at current beta and b
			sqrtrWt = "numeric",# sqrt of weights used with residuals
			sqrtXWt = "numeric",# sqrt of model matrix row weights
			RZX = "dMatrix", # dgeMatrix or dgCMatrix
			RX = "CholKind",  # "Cholesky" (dense) or "dCHMsimpl" (sparse)
			XtX = "dsC_or_dpoMatrix", # "dsymmetricMatrix", # dpo* or dsC*
			Vtr = "numeric"
			),
	 validity = function(object) TRUE ## FIXME .Call(mer_validate, object)
	 )

setClass("merTrms",
         representation(
                        flist = "data.frame", # list of grouping factors
                        cnms = "list",        # of character - column names
                        ncTrms = "integer"    # := sapply(cnms, length)
                        ),
         contains = "mer")
setClass("lmer", representation(),# <-- those that are unique to lmer() case
         contains = "merTrms")

setClass("glmer",
	 representation(# slots	 unique to glmer() :
			eta = "numeric", # unbounded predictor			- GLM
			muEta = "numeric",# d mu/d eta evaluated at current eta - GLM
			var = "numeric"	 # conditional variances of Y		- GLM
			## ghx = "numeric", # zeros of Hermite polynomial
			## ghw = "numeric", # weights used for AGQ
			),
	 contains = "merTrms")

setClass("nlmer",
	 representation(nlmodel = "call" # nonlinear model call
			),
	 contains = "merTrms")


setClass("summary.mer", # Additional slots in a summary object
         representation(
			methTitle = "character",
			logLik= "logLik",
			ngrps = "integer",
			sigma = "numeric", # scale, non-negative number
			coefs = "matrix",
			## vcov = "dpoMatrixorNULL", was  vcov = "dpoMatrix",
			REmat = "matrix",
			AICtab= "data.frame"),
         contains = "mer")
## FIXME -- need differentiate summary.lmer, summary.glmer, .... ?


setClass("glmerenv", contains = "merenvtrms",
         validity = function(object) {
             rho <- env(object)
             n <- length(rho$y)
             if (!(is.numeric(rho$mu) && length(rho$mu) == n))
                 return("environment must contain a numeric vector \"mu\"")
             if (!(is.numeric(rho$muEta) && length(rho$muEta) == n))
                 return("environment must contain a numeric vector \"muEta\"")
             if (!(is.numeric(rho$var) && length(rho$var) == n))
                 return("environment must contain a numeric vector \"var\"")
             if (!(is.numeric(rho$sqrtrWt) && length(rho$sqrtrWt) == n))
                 return("environment must contain a numeric vector \"sqrtrWt\"")
             if (!(is.list(family <- rho$family)))
                 return("environment must contain a list \"family\"")
             if (!(is.character(family$family) && length(family$family) == 1))
                 return("family list must contain a character variable \"family\"")
             if (!(is.character(family$link) && length(family$link) == 1))
                 return("family list must contain a character variable \"link\"")
             TRUE
         })

##' Random-effects module.
##'
##' Zt is the transpose of the sparse model matrix.  The number of
##' rows in Zt may be a multiple of the number of columns in Ut.

setClass("reModule",
         representation(L = "CHMfactor",
                        Lambda = "dgCMatrix",
                        Lind = "integer",
                        Ut = "dgCMatrix",## U := Z Lambda; Ut := U' = Lambda' Z' = Lambda' Zt
                        Zt = "dgCMatrix",## = Z'
                        lower = "numeric",
                        theta = "numeric",
                        u = "numeric",
                        ldL2 = "numeric"),
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
             if (length(object@ldL2) != 1L)
                 return("ldL2 must have length 1")
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
setClass("feModule",
         representation(beta = "numeric", ldRX2 = "numeric", "VIRTUAL"))

##' Dense fixed-effects module
setClass("deFeMod",
         representation(X = "dgeMatrix",
                        RZX = "dgeMatrix",
                        RX = "Cholesky"),
         contains = "feModule",
         validity = function(object) {
             p <- ncol(object@X)
             if (ncol(object@RZX) != p || ncol(object@RX) != p)
                 return("Number of columns in X, RZX and RX must match")
             TRUE
         })

##' Sparse fixed-effects module
setClass("spFeMod",
         representation(X = "dgCMatrix",
                        RZX = "dgCMatrix",
                        RX = "CHMfactor"),
         contains = "feModule",
         validity = function(object) {
             p <- ncol(object@X)
             if (ncol(object@RZX) != p || ncol(object@RX) != p)
                 return("Number of columns in X, RZX and RX must match")
             if (isLDL(object@RX))
                 return("RX must be an LL factor, not LDL")
             TRUE
         })

##' lmer dense fixed-effects module
##'
##' An lmer fixed-effects module is not reweightable so products ZtX
##' and XtX are pre-computed and stored
##'
setClass("lmerDeFeMod",
         representation(ZtX = "dgeMatrix",
                        XtX = "dpoMatrix"),
         contains = "deFeMod",
         validity = function(object) {
             if (any(dim(object@ZtX) != dim(object@RZX)))
                 return("dimensions of ZtX and RZX must match")
             if (any(dim(object@XtX) != dim(object@RX)))
                 return("dimensions of XtX and RX must match")
             if (length(object@ldRX2) != 1L)
                 return("ldRX2 must have length 1")
             TRUE
         })

##' lmer sparse fixed-effects module
##'
setClass("lmerSpFeMod",
         representation(ZtX = "dgCMatrix",
                        XtX = "dsCMatrix"),
         contains = "spFeMod",
         validity = function(object) {
             if (any(dim(object@ZtX) != dim(object@RZX)))
                 return("dimensions of ZtX and RZX must match")
             if (any(dim(object@XtX) != dim(object@RX)))
                 return("dimensions of XtX and RX must match")
             if (length(object@ldRX2) != 1L)
                 return("ldRX2 must have length 1")
             TRUE
         })


##' reweightable dense fixed-effects module
##'
##' a reweightable fixed-effects module contains V, a weighted version
##' of X but no precomputed products.  The number of rows in X can be
##' multiple of the number of rows in V
setClass("rwDeFeMod",
         representation(V = "dgeMatrix"),
         contains = "deFeMod",
         validity = function(object) {
             if (ncol(object@V) != ncol(object@X))
                 return("number of columns in X and V must match")
             TRUE
         })

##' reweightable sparse fixed-effects module
##'
setClass("rwSpFeMod",
         representation(V = "dgCMatrix"),
         contains = "spFeMod",
         validity = function(object) {
             if (ncol(object@V) != ncol(object@X))
                 return("number of columns in X and V must match")
             TRUE
         })

##' mer response module
##' y, offset and mu are as expected
##' wrss is the scalar weighted residual sum of squares
##'
##' Utr is the q-dimensional product of Ut (reModule) and the weighted
##' residuals.
##'
##' Vtr is the p-dimension crossproduct of V (reweightable feModule)
##' or X (lmer feModule) and the weighted residual, which initially is
##' y in the lmer case.
##'
##' cu is the intermediate solution for u, cbeta for beta
setClass("merResp",
         representation(y = "numeric",
                        weights = "numeric", # prior weights
                        offset = "numeric",
                        mu = "numeric",
                        wtres = "numeric", # weighted residuals
                        wrss = "numeric",  # weighted residual sum of squares
                        Utr = "numeric",
                        Vtr = "numeric",
                        cu = "numeric",
                        cbeta = "numeric"),
         validity = function(object) {
             n <- length(object@y)
             if (!(length(object@offset) %in% c(0L, n)))
                 return("length(offset) must be 0 or length(y)")
             if (!(length(object@weights) %in% c(0L, n)))
                 return("length(weights) must be 0 or length(y)")
             if (length(object@mu) != n || length(object@wtres) != n)
                 return("lengths of mu and wtres must match length(y)")
             if (length(object@wrss) != 1L)
                 return("length of wrss must be 1")
             if (length(object@Utr) != length(object@cu))
                 return("lengths of Utr and cu must match")
             if (length(object@Vtr) != length(object@cbeta))
                 return("lengths of Vtr and cbeta must match")
             TRUE
         })

##' reweightable response module
##'
##' gamma is the linear predictor, which is transformed to mu
setClass("rwResp",
         representation(gamma = "numeric",
                        sqrtrwt = "numeric",
                        sqrtXwt = "dgeMatrix"),
         contains = "merResp",
         validity = function(object) {
             lg <- length(object@gamma)
             n <- length(object@y)
             if (lg < 1 || lg %% n)
                 return("length(gamma) must be a positive multiple of length(y)")
             if (length(object@sqrtXwt) != lg)
                 return("length(sqrtXwt) != length(gamma)")
             if (length(object@sqrtrwt) != n)
                 return("length(sqrtrwt) != length(y)")
             TRUE
         })

##' glmer response module
setClass("glmerResp",
         representation(family = "family",
                        muEta = "numeric",
                        n = "numeric",    # for evaluation of the aic
                        var = "numeric"), # variances of responses
         contains = "rwResp",
         validity = function(object) {
             n <- length(object@y)
             lXwt <- length(object@sqrtXwt)
             if (lXwt < 1L || lXwt %% n)
                 return("length(sqrtXwt) must be a positive multiple of length(y)")
             if (length(object@muEta) != n || length(object@var) != n)
                 return("lengths of muEta and var must match length(y)")
         })

##' nlmer response module
setClass("nlmerResp",
         representation(nlenv = "environment",
                        gradient = "matrix"),
         contains = "rwResp",
         validity = function(object) {
             n <- length(object@y)
             N <- length(object@gamma)
             s <- N %/% n
             if (dim(gradient) != c(n, s))
                 return("dimension mismatch on gradient, y and gamma")
         })

##' nglmer response module
setClass("nglmerResp", contains = c("glmerResp", "nlmerResp"))

setClass("lmerMod",
         representation(call = "call",
			frame = "data.frame", # "model.frame" is not S4-ized yet
                        re = "reModule",
                        resp = "merResp",
                        REML = "logical", "VIRTUAL"),
         validity = function(object)
     {
         if (is(object@resp, "rwResp"))
             return("lmer modules cannot be reweightable")
     })

setClass("lmerDe", representation(fe = "lmerDeFeMod"), contains = "lmerMod")

setClass("lmerSp", representation(fe = "lmerSpFeMod"), contains = "lmerMod")

setClass("glmerMod",
         representation(call = "call", frame = "data.frame",
                        re = "reModule",
                        resp = "glmerResp",
                        "VIRTUAL"))

setClass("glmerDe", representation(fe = "rwDeFeMod"), contains = "glmerMod")

setClass("glmerSp", representation(fe = "rwSpFeMod"), contains = "glmerMod")
