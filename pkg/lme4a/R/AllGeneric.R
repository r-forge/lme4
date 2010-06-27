setGeneric("lmList",
           function(formula, data, family, subset, weights,
                    na.action, offset, pool, ...)
           standardGeneric("lmList"))

## utilities, *not* exported (yet ?)
setGeneric("isREML", function(x) standardGeneric("isREML"),
	   valueClass = "logical")
setGeneric("getCall", function(x) standardGeneric("getCall"),
	   valueClass = "call")

## utilities, these *exported*:
setGeneric("getL", function(x) standardGeneric("getL"))

##' extract the deviance components
setGeneric("devcomp", function(x, ...) standardGeneric("devcomp"))


refitML <- function(x) {
    if (!isREML(x)) return(x)
    update(x, REML = FALSE)
}

setGeneric("refitML", function(x) standardGeneric("refitML"))

fixed.effects <- function(object, ...) {
    ## fixed.effects was an alternative name for fixef
    .Deprecated("fixef")
    mCall = match.call()
    mCall[[1]] = as.name("fixef")
    eval(mCall, parent.frame())
}

random.effects <- function(object, ...) {
    ## random.effects was an alternative name for ranef
    .Deprecated("ranef")
    mCall = match.call()
    mCall[[1]] = as.name("ranef")
    eval(mCall, parent.frame())
}

setGeneric("sigma", function(object, ...) standardGeneric("sigma"))

if (FALSE) {
setGeneric("HPDinterval",
           function(object, prob = 0.95, ...) standardGeneric("HPDinterval"))
}

if (FALSE) {
setGeneric("mcmcsamp",
           function(object, n = 1, verbose = FALSE, ...)
           standardGeneric("mcmcsamp"))
}
