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

setClass("optenv", representation(setPars = "function",
                                  getPars = "function",
                                  getBounds = "function"))

setClass("merenv", contains = "optenv",
         validity = function(object)
     {
         rho <- env(object)
         if (!(is.numeric(y <- rho$y) && (n <- length(y)) > 0))
             return("environment must contain a non-trivial numeric response y")
         if (!(is(X <- rho$X, "Matrix") && is(Zt <- rho$Zt, "Matrix") &&
               nrow(X) == ncol(Zt)))
             return("environment must contain Matrix objects X and Zt with nrow(X) == ncol(Zt)")
         if (!(is.numeric(beta <- rho$beta) && length(beta) == ncol(X)))
             return(sprintf("environment must contain a numeric vector beta of length %d", ncol(X)))
         if (!(is.numeric(u <- rho$u) && (q <- length(u)) == nrow(Zt)))
             return(sprintf("environment must contain a numeric vector u of length %d", nrow(Zt)))
         if (!is(Lambda <- rho$Lambda, "Matrix") && all(dim(Lambda) == q))
             return(sprintf("environment must contain a %d by %d Matrix Lambda",
                            q, q))
         TRUE
     })
