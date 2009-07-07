setGeneric("lmList",
           function(formula, data, family, subset, weights,
                    na.action, offset, pool, ...)
           standardGeneric("lmList"))

if (FALSE) {
setGeneric("gsummary",
           function (object, FUN, form, level, groups,
                     omitGroupingFactor = FALSE, 
                     invariantsOnly = FALSE, ...)
           standardGeneric("gsummary"))
}

setGeneric("evalDev", function(object, pars, ...) standardGeneric("evalDev"))

setGeneric("fixef", function(object, ...) standardGeneric("fixef"))

#' returns an A template from Zt
setGeneric("create_A", function(rCF, rho, ...) standardGeneric("create_A"),
           valueClass = "dgCMatrix")

#' returns the parameter bounds
setGeneric("getBounds", function(x, ...) standardGeneric("getBounds"), valueClass = "list")

#' returns Lambda as a square sparse matrix
setGeneric("getLambda", function(rCF, ...) standardGeneric("getLambda"))

#' returns the current parameter values
setGeneric("getPars", function(x, ...) standardGeneric("getPars"), valueClass = "numeric")

#' sets the parameter values
setGeneric("setPars", function(x, pars, ...) standardGeneric("setPars"), valueClass = "numeric")

#' extract the environment associated with an object
setGeneric("env", function(x, ...) standardGeneric("env"), valueClass = "environment")

fixed.effects <- function(object, ...) {
    ## fixed.effects was an alternative name for fixef
    .Deprecated("fixef")
    mCall = match.call()
    mCall[[1]] = as.name("fixef")
    eval(mCall, parent.frame())
}

setGeneric("ranef", function(object, ...) standardGeneric("ranef"))

random.effects <- function(object, ...) {
    ## random.effects was an alternative name for ranef
    .Deprecated("ranef")
    mCall = match.call()
    mCall[[1]] = as.name("ranef")
    eval(mCall, parent.frame())
}

setGeneric("BIC", function(object, ...) standardGeneric("BIC"))

setMethod("BIC", "logLik",
          function(object, ...)
          -2 * (c(object) - attr(object, "df") * log(attr(object, "nobs"))/2)
          )

setGeneric("HPDinterval",
           function(object, prob = 0.95, ...) standardGeneric("HPDinterval"))

setGeneric("mcmcsamp",
           function(object, n = 1, verbose = FALSE, ...)
           standardGeneric("mcmcsamp"))

#setGeneric("pooledSD", function(x, ...) standardGeneric("pooledSD"))

setGeneric("sigma", function(object, ...) standardGeneric("sigma"))

setGeneric("VarCorr", function(x, ...) standardGeneric("VarCorr"))

setGeneric("traceplot", function(x, ...) standardGeneric("traceplot"))

setGeneric("refit", function(object, newresp, ...) standardGeneric("refit"))
