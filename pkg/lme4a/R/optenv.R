#' Optimization-related methods for the environment class

setMethod("env", "optenv", function(x, ...) environment(x@getPars))

## In newer versions of R,  the method argument list is *extended*
## to the one of the generic
if(getRversion() >= "2.11.0") {
setMethod("nlminb", signature(start = "optenv"),
          function(start, ...)
      {
          ## the argument name 'control' is matched from the generic
          control <- as.list(control)
          bb <- start@getBounds()
          nlminb(start@getPars(), start@setPars, lower = bb[,"lower"],
                 upper = bb[, "upper"], control = control)
      })

setMethod("bobyqa", signature(par = "optenv"),
          function(par, ...)
      {
          ## the argument name 'control' is matched from the generic
          control <- as.list(control)
          bb <- par@getBounds()
          bobyqa(par@getPars(), par@setPars, lower = bb[,"lower"],
                 upper = bb[, "upper"], control = control)
      })

setMethod("optimize", signature(f = "optenv"),
          function(f, ...)
      {
          ## the argument name 'tol' is matched from the generic
          stopifnot(length(f@getPars()) == 1)
          bb <- f@getBounds()
          ## hmm, typically have  (0, Inf) here,  but optimize really starts at
          ##  a + (1-phi)*(b-a) = phi*a + (1-phi)b  = .382*a + .618*b,
          ## and its extremely silly to have  (a = 0, b = 1.8e308)
	  if(any(is.infinite(bb))) ## optimize does not allow '+- Inf' bounds:
	      bb[] <- pmin(.Machine$double.xmax,
			   pmax(-.Machine$double.xmax, bb))
          r <- optimize(f@setPars, lower = bb[,"lower"], upper = bb[,"upper"], tol=tol)
      })

setMethod("optimize", signature(f = "merenv"),
          function(f, ...)
      {
          ## the argument name 'tol' is matched from the generic
          stopifnot(length(f@getPars()) == 1)
          optimize(f@setPars, lower = 0, upper = 100, tol=tol)
      })

} else {## R <= 2.10.x -- need correct (full) argument list in method

setMethod("nlminb", signature(start = "optenv"),
          function (start, objective, gradient = NULL, hessian = NULL, ...,
                    scale = 1, control = list(), lower = -Inf, upper = Inf)
      {
          control <- as.list(control)
          bb <- start@getBounds()
          nlminb(start@getPars(), start@setPars, lower = bb[,"lower"],
                 upper = bb[, "upper"], control = control)
      })

setMethod("bobyqa", signature(par = "optenv"),
          function(par, fn, lower = -Inf, upper = Inf, control = list(), ...)
      {
          control <- as.list(control)
          bb <- par@getBounds()
          bobyqa(par@getPars(), par@setPars, lower = bb[,"lower"],
                 upper = bb[, "upper"], control = control)
      })

setMethod("optimize", signature(f = "optenv"),
          function (f, interval, ..., lower = min(interval), upper = max(interval),
                    maximum = FALSE, tol = .Machine$double.eps^0.25)
      {
          stopifnot(length(f@getPars()) == 1)
          bb <- f@getBounds()
          optimize(f@setPars, lower = bb[,"lower"], upper = bb[,"upper"], tol=tol)
      })

}## if( R >= 2.11.0 ) .... else ....
