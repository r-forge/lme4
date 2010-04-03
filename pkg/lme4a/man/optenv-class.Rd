\name{optenv-class}
\Rdversion{1.1}
\docType{class}
\alias{optenv-class}% --> ./merenv-class.Rd now has those subclasses
\alias{bobyqa,optenv-method}
\alias{nlminb,optenv-method}
\alias{optimize,optenv-method}
\title{Class "optenv" - an optimization environment}
\description{
  The \code{\link{class}} \code{"optenv"} is an S4 class of functions
  used for optimization of an objective function.  The functions, which
  are the slots of the class, share a common environment in which the
  state of the system being optimized is stored.  The S4 class \code{"lmerenv"}
  inherits from "optenv" and has objects in the environment representing
  a linear mixed-effects model.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("optenv", ...)}.
  Generally it is objects of subclasses, such as
  \code{\linkS4class{merenv}},
  that are created directly.
}
\section{Slots}{
  \describe{
    \item{\code{setPars}:}{A function that sets new parameter values in
      the object and returns the new value of the objective.}
    \item{\code{getPars}:}{A function that returns the current value of
      the parameters are a numeric vector.}
    \item{\code{getBounds}:}{A function that returns the bounds on the
      parameter values as a numeric matrix with 2 columns,
      \code{"lower"} and \code{"upper"}.  The number of rows must be the
      length of the value of the \code{getPars} slot.}
  }
}
\section{Methods}{
  \describe{
    \item{bobyqa}{\code{signature(par = "optenv")}: optimize the
      objective using \code{bobyqa} from the \pkg{minqa} package.}
    \item{nlminb}{\code{signature(start = "optenv")}: optimize the
      objective using \code{\link{nlminb}}.}
    \item{optimize}{\code{signature(f = "optenv")}: optimize the
      objective using \code{\link{optimize}}. This method only applies
      to functions of a scalar argument.}
  }
}
\seealso{
  \code{\link{nlminb}}, \code{\link{optimize}};
  subclasses of \code{"optenv"}: \code{\linkS4class{merenv}}, etc.
}
\examples{
showClass("optenv")
%% Consider mentioning that different subclasses have different validity
%% methods :
%%  str( getValidity(getClass("merenv")) )

fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
            doFit = FALSE) #- so we get the "merenv" (before optimization)
nlminb(fm1) # optimze it --> *does* change the 'fm1' object !!!!!
fm1.s <- copyMerenv(fm1)# save it
str(n.fm1 <- nlminb(fm1, control= list(trace=TRUE, rel.tol = 1e-11)))
## but did not really get better:
all.equal(fm1@getPars(), n.fm1$par, tol=0)
   # -> *no* difference on [64b AMD Linux 2.6.30]

## copy back:
fm1 <- copyMerenv(fm1.s)

stopifnot(require("minqa")) ## works as we depend on it ..

str(b1.fm1 <- bobyqa(fm1))
b2.fm1 <- bobyqa(fm1, control= list(rhoend = 1e-13, iprint = 1))

stopifnot(  all.equal(b1.fm1[["par"]],  b2.fm1[["par"]],  tol = 1e-6)
          , all.equal(b1.fm1[["fval"]], b2.fm1[["fval"]], tol = 1e-15)
          , b1.fm1$feval < b2.fm1$feval
          , all.equal(b1.fm1$fval, fm1@setPars(b1.fm1$par), tol= 1e-15)
          , all.equal(b2.fm1$fval, fm1@setPars(b2.fm1$par), tol= 1e-15)
          , TRUE)
}
\keyword{classes}
