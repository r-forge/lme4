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
      objective using \code{bobyqa} from the \code{minqa} package.}
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
}
\keyword{classes}
