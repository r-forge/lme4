\name{merenv-class}
\title{Class "merenv" - environments from lmer() etc}
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
\alias{model.frame,merenv-method}
\alias{model.matrix,merenv-method}
\alias{terms,merenv-method}
\alias{update,merenv-method}
\alias{vcov,merenv-method}
\alias{VarCorr,merenv-method}
%
\alias{coef,merenvtrms-method}
% *lmerenv methods: typically lmerenv and merenv (= non-lmerenv)
\alias{anova,lmerenv-method}
\alias{deviance,merenv-method}
\alias{deviance,lmerenv-method}
\alias{logLik,lmerenv-method}
\alias{logLik,merenv-method}
\alias{print,merenv-method}
\alias{print,lmerenv-method}
\alias{show,merenv-method}
\alias{show,lmerenv-method}
\alias{summary,lmerenv-method}
\alias{summary,glmerenv-method}% for now
%
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
  However, typically, these are are created directly, from
  \code{\link{lmer1}()},  \code{\link{glmer1}()}, etc.
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
    \item{model.frame}{\code{signature(formula = "merenv")}: ... }
    \item{model.matrix}{\code{signature(object = "merenv")}: ... }
    \item{sigma}{\code{signature(object = "merenv")}: ... }
    \item{terms}{\code{signature(x = "merenv")}: ... }
    \item{update}{\code{signature(object = "merenv")}: ... }
    \item{VarCorr}{\code{signature(object = "merenv")}: Returns a
  summary of the variance-component parameters. }
    \item{vcov}{\code{signature(object = "merenv")}: Returns the
  variance-covariance matrix of the fixed-effects parameters. }

    %% merenvtrms :
    \item{coef}{\code{signature(object = "merenvtrms")}: ... }

    %% lmerenv :
    \item{anova}{\code{signature(object = "lmerenv")}: ... }
    \item{deviance}{\code{signature(object = "lmerenv")}: ... }
    \item{logLik}{\code{signature(object = "lmerenv")}: ... }
    \item{print}{\code{signature(x = "lmerenv")}: ... }
    \item{show}{\code{signature(object = "lmerenv")}: ... }
    \item{summary}{\code{signature(object = "lmerenv")} and additional
      argument \code{varcov = FALSE}:  The resulting (S3 class)
      \code{"summary.lmer"} object is \code{\link{list}}-like, with a
      non-\code{NULL} component \code{vcov} only when explicitly requested
      by \code{summary(..., varcov = TRUE)}.}
  }
}
% \references{ ~put references to the literature/web site here ~ }
\seealso{
  \code{\link{env}()} as a shortcut access to the environment with all the
  objects.
  the \code{\link{VarCorr}} method.
}
\examples{
showClass("merenv") # note its subclasses
%% Consider mentioning that different subclasses have different validity
%% methods :
%%  str( getValidity(getClass("merenv")) )
}
\keyword{classes}
