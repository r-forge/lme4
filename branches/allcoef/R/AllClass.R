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

setClass("mer", # mixed-effects representation used for PIRLS
	 representation(## original data
                        env = "environment", # evaluation env for nonlinear model
                        nlmodel = "call",    # nonlinear model call
### FIXME: Move frame up in class hierarchy                        
                        frame = "data.frame",# model frame (or empty frame)
### FIXME: Just need fl1 as a integer vector instead of flist?
                        flist = "data.frame",# list of grouping factors
                        X = "matrix",        # fixed effects model matrix
### FIXME: Probably do without this for the time being?
                        Xst = "dgCMatrix",   # sparse fixed effects model matrix
                        Zt = "dgCMatrix",    # sparse form of Z'
                        pWt = "numeric",     # prior weights,
                        offset = "numeric",  # length 0 -> no offset
                        y = "numeric",       # response vector
                        dims = "integer",    # dimensions and indicators
                        ## slots that vary during optimization
                        etaGamma = "matrix", # gradient matrix
                        perm = "integer",    # permutation vector
                        A = "dgCMatrix",     # (Z %*% Lambda)'
                        L = "CHMfactor",     # Cholesky factor of U'U + I
                        deviance = "numeric",# deviance and components
			fixef = "numeric",   # fixed effects (length p)
			ranef = "numeric",   # random effects (length q)
                        u = "numeric",       # orthogonal random effects (q)
                        eta = "numeric",     # unbounded predictor
                        mu = "numeric",      # conditional mean
                        muEta = "numeric",   # d mu/d eta at current eta
                        var = "numeric",     # conditional variances of Y
                        resid = "numeric",   # raw residuals
                        sqrtrWt = "numeric", # sqrt of weights used with residuals
                        RZX = "matrix",      # dense sol. to L RZX = UV
                        RX = "matrix",       # Cholesky factor of downdated V'V
		        ghx = "numeric",     # zeros of Hermite polynomial
			ghw = "numeric"),    # weights used for AGQ
         validity = function(object) .Call(mer_validate, object))

setClass("ST",  # Scale/unit Triangular representation of relative covariance
         representation(ST = "list", # list of TSST' rep of rel. cov. mats
                        Gp = "integer"), # pointers to r.e. term groups
         validity = function(object) .Call(ST_validate, object))

##' Parameterized components of an mer model.
##'
##' The virtual class of parameterized components of an mer model
setClass("merParam", representation("VIRTUAL"))

##' List of merParam objects
setClass("merParamList",
         representation(offset = "integer"), # pointers into parameter vector
         contains = "list",
         validity = function(object)
         all(unlist(lapply(object, "is", class2 = "merParam"))))

##' Random-effects covariance factors.
##'
##' The virtual class of components that generate the factors of the
##' relative covariance matrices for random effects.
setClass("merREfac",
         representation(offset = "integer", # offset into the parameter vector
                        "VIRTUAL"),
         contains = "merParam")

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
