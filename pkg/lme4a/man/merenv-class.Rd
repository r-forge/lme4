\name{merenv-class}
\Rdversion{1.1}
\docType{class}
\alias{merenv-class}
\alias{merenvtrms-class}
\alias{lmerenv-class}
\alias{glmerenv-class}
%
%\alias{fixef,merenv-method} -> ./fixef.Rd
%\alias{ranef,merenv-method} -> ./ranef.Rd
\alias{formula,merenv-method}
\alias{update,merenv-method}
% subclasses methods
\alias{anova,lmerenv-method}
\alias{coef,merenvtrms-method}
\alias{deviance,lmerenv-method}
\alias{logLik,lmerenv-method}
\alias{model.frame,lmerenv-method}
\alias{model.matrix,lmerenv-method}
\alias{print,lmerenv-method}
\alias{profile,lmerenv-method}
\alias{show,lmerenv-method}
\alias{summary,lmerenv-method}
\alias{terms,lmerenv-method}
\alias{VarCorr,lmerenv-method}
\alias{vcov,lmerenv-method}

\title{Class "merenv" - environments from lmer() etc}
\description{
  The \code{\link{class}} \code{"merenv"} is extending (aka
  \sQuote{inheriting from}) \code{"\linkS4class{optenv}"}, an S4 class
  of functions.  The functions, which are the slots of the class, share
  a common environment representing the state of the system being
  optimized. Apart from the optimization state, it has objects in the
  environment representing \bold{m}ixed-\bold{e}ffects models in \R,
  hence \dQuote{mer}.
}
\section{Objects from the Class}{
  \code{"merenv"} is a virtual class.  Objects of (\dQuote{actual})
  subclasses can be created however, by calls of the form
  \code{new("lmerenv", ...)}.
  However, typically, these are are created directly, e.g., from
  \code{\link{lmer}()}.
}
\section{Slots}{
  All slots are inherited from \code{"\linkS4class{optenv}"}.
}
\section{Extends}{
  Class \code{"merenv"} extends class \code{"\linkS4class{optenv}"},
  directly.

  Class \code{"lmerenv"} extends class \code{"merenvtrms"}, directly,
  and
  classes \code{"merenv"}, by class \code{"merenvtrms"}, distance 2;
  and class \code{"\linkS4class{optenv}"}, by class \code{"merenvtrms"},
  distance 3.
}
\section{Methods}{
  \describe{
    %% merenv :
    \item{formula}{\code{signature(x = "merenv")}: ... }
    \item{update}{\code{signature(object = "merenv")}: ... }

    %% merenvtrms :
    \item{coef}{\code{signature(object = "merenvtrms")}: ... }

    %% lmerenv :
    \item{anova}{\code{signature(object = "lmerenv")}: ... }
    \item{deviance}{\code{signature(object = "lmerenv")}: ... }
    \item{logLik}{\code{signature(object = "lmerenv")}: ... }
    \item{model.frame}{\code{signature(formula = "lmerenv")}: ... }
    \item{model.matrix}{\code{signature(object = "lmerenv")}: ... }
    \item{print}{\code{signature(x = "lmerenv")}: ... }
    \item{profile}{\code{signature(fitted = "lmerenv")}: ... }
    \item{show}{\code{signature(object = "lmerenv")}: ... }
    \item{sigma}{\code{signature(object = "lmerenv")}: ... }
    \item{summary}{\code{signature(object = "lmerenv")}: ... }
    \item{terms}{\code{signature(x = "lmerenv")}: ... }
    \item{VarCorr}{\code{signature(x = "lmerenv")}: ... }
    \item{vcov}{\code{signature(object = "lmerenv")}: ... }
  }
}
% \references{ ~put references to the literature/web site here ~ }
% \author{ ~~who you are~~ }
\seealso{
  \code{\link{env}()} as a shortcut access to the environment with all the
  objects.
  \code{\link{getPars}}, \code{\link{getBounds}} and \code{\link{setPars}}.
}
\examples{
showClass("merenv") # note its subclasses
%% Consider mentioning that different subclasses have different validity
%% methods :
%%  str( getValidity(getClass("merenv")) )
}
\keyword{classes}
