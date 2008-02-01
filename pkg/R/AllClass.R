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

setClass("mer",
	 representation(## original data
#                        famName = "character",# name of GLM family and link
                        env = "environment",# evaluation env for nonlinear model
                        nlmodel = "call",# nonlinear model call
                        frame = "data.frame",# model frame (or empty frame)
                        call = "call",   # matched call
                        flist = "data.frame",  # list of grouping factors
                        X = "matrix",    # fixed effects model matrix
                        Zt = "dgCMatrix",# sparse form of Z'
                        pWt = "numeric",# prior weights,
                        offset = "numeric", # length 0 -> no offset
                        y = "numeric",   # response vector
                        cnames = "list", # row/column names of els of ST
                        Gp = "integer",  # pointers to row groups of Zt
                        dims = "integer",# dimensions and indicators
                        ## slots that vary during optimization
                        ST = "list", # list of TSST' rep of rel. var. mats
                        V = "matrix",    # gradient matrix
                        A = "dgCMatrix", # (ZTS)'
                        Cm = "dgCMatrix", # AH'G^{-1}W^{1/2} when s > 0
                        Cx = "numeric",  # x slot of Cm when s == 1 (full Cm not stored)
                        L = "CHMfactor", # Cholesky factor of weighted P(AA' + I)P'
                        deviance = "numeric", # ML and REML deviance and components
			fixef = "numeric",# fixed effects (length p)
			ranef = "numeric",# random effects (length q)
                        u = "numeric",   # orthogonal random effects (q)
                        eta = "numeric", # unbounded predictor
                        mu = "numeric",  # fitted values at current beta and b
                        muEta = "numeric",# d mu/d eta evaluated at current eta
                        var = "numeric", # conditional variances of Y
                        resid = "numeric",# raw residuals at current beta and b
                        sqrtXWt = "matrix",# sqrt of model matrix row weights
                        sqrtrWt = "numeric",# sqrt of weights used with residuals
                        RCXy = "matrix", # dense sol. to L RCXy = S T'ZtXy
                        RXy = "matrix"), # Cholesky factor of downdated XytXy
         validity = function(object) .Call(mer_validate, object))

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

setClass("ranef.mer", contains = "list")

setClass("coef.mer", contains = "list")

setClass("pedigree", representation =
	 list(sire = "integer", dam = "integer", label = "character"),
	 validity = function(object) {
	     n <- length(sire <- object@sire)
	     if (length(dam <- object@dam) != n)
		 return("sire and dam slots must be the same length")
	     if (length(object@label) != n)
		 return("'label' slot must have the same length as 'sire' and 'dam'")
	     if(n == 0) return(TRUE)
	     animal <- 1:n
	     snmiss <- !is.na(sire)
	     dnmiss <- !is.na(dam)
	     if (any(sire[snmiss] >= animal[snmiss]) ||
		 any(dam[dnmiss] >= animal[dnmiss]))
		 return("the sire and dam must precede the offspring")
             if (any(sire[snmiss] < 1 | sire[snmiss] > n) |
                 any(dam[dnmiss] < 1 | dam[dnmiss] > n))
                 return(paste("Non-missing sire or dam must be in [1,",
                              n, "]", sep = ''))
	     TRUE
	 })
