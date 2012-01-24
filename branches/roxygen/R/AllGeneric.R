
## utilities, these *exported*:
##' @export getL
setGeneric("getL", function(x) standardGeneric("getL"))

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

##' Extract the residual standard error from a fitted model.
##'
##' This is a generic function.  At present the only methods are for mixed-effects
##' models of class \code{\linkS4class{merMod}}.
##' @title Extract residual standard error
##' @param object a fitted model.
##' @param ... additional, optional arguments.  (None are used in the merMod method)
##' @return the residual standard error as a scalar
##' @export
sigma <- function(object, ...) UseMethod("sigma")

##' Check if a model has been fit according to the REML criterion
##'
##' This is a generic function.  At present the only methods are for mixed-effects
##' models of class \code{\linkS4class{merMod}}.
##' @title Check for a REML fit
##' @param x a fitted model.
##' @param ... additional, optional arguments.  (None are used in the merMod method)
##' @return \code{TRUE} if \code{x} has been fit by REML, otherwise FALSE
##' @export
isREML <- function(x, ...) UseMethod("isREML")

##' Refit a model using the maximum likelihood criterion
##'
##' This function is primarily used to get a maximum likelihood fit of
##' a linear mixed-effects model for an \code{\link{anova}} comparison.
##' @title Refit a model by maximum likelihood criterion
##' @param x a fitted model, usually of class \code{"\linkS4class{lmerMod}"},
##'     to be refit according to the maximum likelihood criterion
##' @param ... optional additional parameters.  None are used at present.
##' @return an object like \code{x} but fit by maximum likelihood
##' @export
refitML <- function(x, ...) UseMethod("refitML")

if (FALSE) {
setGeneric("HPDinterval",
           function(object, prob = 0.95, ...) standardGeneric("HPDinterval"))
}

if (FALSE) {
setGeneric("mcmcsamp",
           function(object, n = 1, verbose = FALSE, ...)
           standardGeneric("mcmcsamp"))
}
