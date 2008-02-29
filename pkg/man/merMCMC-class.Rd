\name{merMCMC-class}
\docType{class}
% Classes
\alias{merMCMC-class}
% Methods
%\alias{anova,mer-method}
\title{Mixed-model Markov chain Monte Carlo results}
\description{Elements of the \code{merMCMC} class are Markov chain Monte
  Carlo samples from a fitted mixed-effects model.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("merMCMC", ...)}
  or, more commonly via the \code{mer} method for the generic
  \code{\link{mcmcsamp}} function.
}
\section{Slots}{
  \describe{
    \item{\code{Gp}:}{a copy of the \code{Gp} slot of the original
      \code{\linkS4class{mer}} object}
    \item{\code{ST}:}{matrix of samples from the parameters determining
      the \code{ST} slot of the \code{\linkS4class{mer}} object}
    \item{\code{call}:}{Matched call from the original
      \code{\linkS4class{mer}} object}
    \item{\code{deviance}:}{vector of samples of the (ML) deviance}
    \item{\code{dims}:}{a copy of the \code{dims} slot of the original
      \code{\linkS4class{mer}} object}
    \item{\code{fixef}:}{matrix of samples of the fixed-effects parameters}
    \item{\code{nc}:}{Integer vector of length \code{dims["nf"]}.  The
      number of columns of random effects in each term.}
    \item{\code{ranef}:}{matrix of samples of the random effects.  This
      matrix has zero columns unless \code{saveb = TRUE} is specified in
      the call to \code{\link{mcmcsamp}}.  Consider the size of this
      matrix, which could be very large, before setting \code{saveb = TRUE}.}
    \item{\code{sigma}:}{vector of samples of the common scale
      parameter or \code{numeric(0)} if \code{dims["useSc"]} is \code{FALSE}.}
  }
}
\section{Methods}{
  No methods defined with class "merMCMC" in the signature.
}
\seealso{
  \code{\link{mcmcsamp}} produces these objects,
  \code{\link{lmer}}, \code{\link{glmer}} and \code{\link{nlmer}}
  produce the \code{\linkS4class{mer}} objects.
}
\examples{
showClass("merMCMC")
}
\keyword{classes}
