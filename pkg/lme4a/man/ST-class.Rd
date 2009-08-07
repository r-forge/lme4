\name{ST-class}
\Rdversion{1.1}
\docType{class}
\alias{ST-class}
\alias{chol,ST-method}
%\alias{getBounds,ST-method}
%\alias{getPars,ST-method}
\alias{ranef,ST-method}
%\alias{setPars,ST,ANY-method}

\title{Class "ST" - relative covariance factors}
\description{The relative covariance matrix for the random effects in a
  mixed-effects model is a symmetric, positive semidefinite matrix of
  size \eqn{q\times q}{q x q} where \eqn{q} is the number of random
  effects.  Because it is positive semidefinite it can be expressed as a
  factor, which we write as \eqn{\Lambda}, multiplied by the transpose
  of \eqn{\Lambda}.  The factor, \eqn{\Lambda}, is called the relative
  covariance factor.  In the \code{"ST"} class the factor is expressed
  as the product of a unit lower triangular matrix, \eqn{T}, and a
  diagonal matrix, \eqn{S}, with non-negative diagonal elements.}
\section{Objects from the Class}{
  Although objects can be created by calls of the form \code{new("ST",
    ...)} they are typically created in a model-fitting function such as
  \code{\link{lmer}}, \code{\link{glmer}} or \code{\link{nlmer}}.
}
\section{Slots}{
  \describe{
    \item{\code{ST}:}{Object of class \code{"list"} }
    \item{\code{Gp}:}{Object of class \code{"integer"}
      representing the group pointers.  (This will change in
      future releases.) }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{reCovFac}"}, directly.
}
\section{Methods}{
  \describe{
    \item{chol}{\code{signature(x = "ST")}: return the sparse Cholesky
      factor, \eqn{\Lambda}}
    \item{ranef}{\code{signature(object = "ST")}: return the random
      effects in the original scale.}
  }
}
%\references{}
%\author{}
%\note{}
\seealso{\code{\link{lmer}}, \code{\link{glmer}}, \code{\link{nlmer}}
  and \code{\linkS4class{mer}}.
}
\examples{
showClass("ST")
}
\keyword{classes}
