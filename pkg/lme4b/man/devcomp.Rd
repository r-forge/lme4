\name{devcomp}
\title{Extract the components of the deviance}
\alias{devcomp}
\alias{devcomp-methods}
\alias{devcomp,glmerenv-method}
\alias{devcomp,merMod-method}
\alias{devcomp,lmerenv-method}
\alias{devcomp,mer-method}
\description{
  Extract the components of the deviance from a \code{*lmer()} fit.
  This is a generic function with methods for all kinds of classes
  returned from calls to \code{\link{lmer1}}, \code{\link{glmer1}}, etc.
}
\usage{
devcomp(x, ...)
}
\arguments{
  \item{x}{a model fit by \code{\link{lmer1}}.}
  \item{\dots}{optional additional arguments.  None are used at present.}
}
\details{
  Extract several numeric quantities used to evaluate the deviance or
  the REML criterion of the model.
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "glmerenv")}}{ }%% FIXME
    \item{\code{signature(x = "glmerMod")}}{ }

    \item{\code{signature(x = "lmerenv")}}{ }
    \item{\code{signature(x = "lmerMod")}}{ }

    \item{\code{signature(x = "mer")}}{ }
}}
\value{a list of two named vectors
  \item{cmp }{Numeric components including \code{ldL2}, the logarithm of
    the square of the determinant of the sparse Cholesky factor,
    \code{L}; \code{ldRX2}, the logarithm of the square of the
    determinant of the sparse or dense Cholesky factor \eqn{{\boldmath
	R}_X}{RX}; \code{prss}, the penalized residual sum of squares;
    \code{deviance}, the deviance; and \code{REML}, the REML criterion}
  \item{dims}{An integer vector of dimensions including \code{n}, the
    dimension of the response vector; \code{p}, the dimension of the
    fixed-effects parameter vector; \code{nmp}, \eqn{n-p}; and \code{q},
    the dimension of the random effects vector.}
}
\seealso{\code{\link{lmer1}}, \code{\linkS4class{lmerenv}}}
\examples{
devcomp(lmer1(Yield ~ 1 + (1|Batch), Dyestuff))
}
\keyword{models}
\keyword{methods}
