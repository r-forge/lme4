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
         contains = "reCovFac",
         validity = function(object) .Call(ST_validate, object))

setClass("PIRLS", # penalized iteratively reweighted least squares
	 representation(## original data
                        rho = "environment", # evaluation environment for PIRLS
                        nlenv = "environment", # evaluation env for nonlinear model
                        nlmodel = "call",    # nonlinear model call
                        X = "matrix",        # fixed effects model matrix
                        pWt = "numeric",     # prior weights,
                        offset = "numeric",  # length 0 -> no offset
                        y = "numeric",       # response vector
                        dims = "integer"),    # dimensions and indicators
#                        fl1 = "integer",     # first grouping factor
#                        ghw = "numeric",     # weights used for AGQ
#                        ghx = "numeric"),     # zeros of Hermite polynomial
                        ## slots that vary during optimization
#                        A = "dgCMatrix",     # (Z %*% Lambda)'
#                        L = "CHMfactor",     # Cholesky factor of U'U + I
#                        RX = "matrix",       # Cholesky factor of downdated V'V
#                        RZX = "matrix",      # dense sol. to L RZX = UV
#                        deviance = "numeric",# deviance and components
#                        eta = "numeric",     # unbounded predictor
#                        etaGamma = "matrix", # gradient matrix
#                        fixef = "numeric",   # fixed effects (length p)
#                        mu = "numeric",      # conditional mean
#                        muEta = "numeric",   # d mu/d eta at current eta
#                        resid = "numeric",   # raw residuals
#                        sqrtrWt = "numeric", # sqrt of weights used with residuals
#                        u = "numeric",       # orthogonal random effects (q)
#                        var = "numeric"),    # conditional variances of Y
         validity = function(object) .Call(mer_validate, object))

setClass("mer",
         representation(rCF = "reCovFac",
                        PLS = "PIRLS",
#                        Zt = "dgCMatrix",    # sparse form of Z'
#                        flist = "data.frame",# list of grouping factors
                        frame = "data.frame",# model frame
                        ranef = "numeric")   # random effects (length q)
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
