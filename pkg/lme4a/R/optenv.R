#' Optimization-related methods for the environment class

if (FALSE) {
setMethod("getPars", "optenv", function(x, ...) x@getPars())

setMethod("getBounds", "optenv", function(x, ...) x@getBounds())

setMethod("setPars", representation(x = "optenv", pars = "numeric"), function(x, pars, ...) x@setPars(pars))
}

setMethod("env", "optenv", function(x, ...) environment(x@getPars))

setMethod("nlminb", signature(start = "optenv"),
          function(start, ...)
      {
          control <- as.list(list(...)$control)
          bb <- start@getBounds()
          nlminb(start@getPars(), start@setPars, lower = bb[,"lower"],
                 upper = bb[, "upper"], control = control)
      })

setMethod("bobyqa", signature(par = "optenv"),
          function(par, ...)
      {
          control <- do.call(bobyqa.control, as.list(list(...)$control))
          bb <- par@getBounds()
          bobyqa(par@getPars(), par@setPars, lower = bb[,"lower"],
                 upper = bb[, "upper"])
      })

setMethod("optimize", signature(f = "optenv"),
          function(f, ...)
      {
          stopifnot(length(f@getPars()) == 1)
          bb <- f@getBounds()
          optimize(f@setPars, lower = bb[,"lower"], upper = bb[,"upper"])
      })
