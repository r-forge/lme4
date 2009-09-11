## Class definitions for the package

setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

setClass("lmList.confint", contains = "array")

## -------------------- lmer-related Classes --------------------------------

setOldClass("data.frame")
setOldClass("family")
setOldClass("logLik")

setClass("reCovFac",                    # random-effects covariance factor
         representation("VIRTUAL"))

setClass("ST",  # Scale/unit Triangular representation of relative covariance
         representation(ST = "list", # list of TSST' rep of rel. cov. mats
                        Gp = "integer"), # pointers to r.e. term groups
         contains = "reCovFac"
         ## , validity = function(object) .Call(ST_validate, object)
         )

setClass("mer",
	 representation(                       # original data
                        nlenv = "environment",
                        nlmodel = "call", # nonlinear model call
                        frame = "data.frame", # model frame
                        call = "call",        # matched call
                        flist = "data.frame", # list of grouping factors
                        X = "matrix",    # fixed effects model matrix
                        Zt = "dgCMatrix",# sparse form of Z'
                        pWt = "numeric",# prior weights,
                        offset = "numeric", # length 0 -> no offset
                        y = "numeric",   # response vector
                        dims = "integer",# dimensions and indicators
                        ## derived slots
                        rCF = "reCovFac", # relative covariance factor
                        etaGamma = "matrix", # gradient matrix
                        A = "dgCMatrix",     # (Z Gamma)'
                        L = "CHMfactor", # Cholesky factor of U'U+I
                        perm = "integer", # permutation of u
                        deviance = "numeric", # ML and REML deviance and components
			fixef = "numeric",# fixed effects (length p)
			ranef = "list",   # random effects (list of matrices)
                        u = "numeric", # orthogonal random effects (q)
                        eta = "numeric", # unbounded predictor
                        mu = "numeric",  # fitted values at current beta and b
                        muEta = "numeric",# d mu/d eta evaluated at current eta
                        var = "numeric", # conditional variances of Y
                        resid = "numeric",# raw residuals at current beta and b
                        sqrtrWt = "numeric",# sqrt of weights used with residuals
                        RZX = "matrix", # dense sol. to L RZX = ST'ZtX = AX
                        RX = "matrix",  # Cholesky factor of downdated X'X
		        ghx = "numeric", # zeros of Hermite polynomial
			ghw = "numeric") # weights used for AGQ
#         ,validity = function(object) .Call(mer_validate, object)
         )

setClass("merMCMC",
         representation(
                        Gp = "integer",   # Gp slot from mer object
                        ST = "matrix",    # matrix of sampled ST pars
                        call = "call",    # matched call
                        deviance = "numeric",# vector of sampled deviances
                        dims = "integer", # dims from original mer object
                        fixef = "matrix", # matrix of sampled fixed effects pars
                        nc = "integer",   # number of columns per r.e. term
                        ranef = "matrix", # optional matrix of sampled r.e.
                        sigma = "matrix"  # sigma samples (may have 0 columns)
                        ),
         validity = function(object) .Call(merMCMC_validate, object))

setClass("summary.mer",                 # Additional slots in a summary object
         representation(
			methTitle = "character",
			logLik= "logLik",
			ngrps = "integer",
			sigma = "numeric", # scale, non-negative number
			coefs = "matrix",
			vcov = "dpoMatrix",
			REmat = "matrix",
			AICtab= "data.frame"),
         contains = "mer")

#setClass("ranef.mer", contains = "list")

#setClass("coef.mer", contains = "list")

setClass("sparseRasch", representation =
         list(dims = "integer",
              Zt = "dgCMatrix",
              y = "numeric",
              deviance = "numeric",
              offset = "numeric",
              L = "CHMfactor",
              fixef = "numeric",
              mu = "numeric",
              muEta = "numeric",
              pWt = "numeric",
              resid = "numeric",
              sqrtrWt = "numeric",
              var = "numeric"),
         validity = function(object) TRUE)

setClass("merExt",
         representation(X0 = "matrix",    # original fixed effects model matrix
                        Zt0 = "dgCMatrix",# original sparse form of Z'
                        pars = "numeric", # additional parameters
                        y0 = "numeric"),  # original response vector
         contains = "mer"
         )

setClass("lmerStratVar",
         representation(sfac = "factor"),
         contains = "merExt")

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
         validity = function(object)
     {
         rho <- env(object)
         if (!(is.list(flist <- rho$flist) &&
               all(sapply(flist, is.factor))))
             return("environment must contain a list of factors, flist")
         flseq <- seq_along(flist)
         if (!(is.integer(asgn <- attr(flist, "assign")) &&
               all(flseq %in% asgn) &&
               all(asgn %in% flseq)))
             return("asgn attribute of flist missing or malformed")
         nl <- sapply(flist, function(x) length(levels(x)))[asgn]
         if (!(is.list(cnms <- rho$cnms) &&
               all(sapply(cnms, is.character)) &&
               all(sapply(cnms, length) > 0) &&
               length(cnms) == length(asgn)))
             return("list of column names, cnms, must match asgn attribute in length")
     })

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
         if (!(is(RX <- rho$RX, "dMatrix") &&
               all(dim(RX) == c(p, p))))
             return("environment must contain a p by p Matrix RX")
         if (!(is.numeric(Xty <- rho$Xty) &&
               length(Xty == p)))
             return("environment must contain a numeric p-vector Xty")
         if (!(is(XtX <- rho$XtX, "dMatrix") &&
               is(XtX, "symmetricMatrix") &&
               all(dim(RZX) == c(q, p))))
             return("environment must contain a symmetric Matrix XtX")
         TRUE
     })

