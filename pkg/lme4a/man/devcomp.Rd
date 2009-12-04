\name{devcomp}
\alias{devcomp}
\title{Extract the components of the deviance}
\description{
  Extract the components of the deviance from an
  \code{\linkS4class{lmerenv}} object.
}
\usage{
devcomp(x, ...)
}
\arguments{
  \item{x}{a model fit by \code{\link{lmer}}.}
  \item{\dots}{optional additional arguments.  None are used at present.}
}
\details{
  Extract several numeric quantities used to evaluate the deviance or
  the REML criterion of the model.
}
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
\seealso{\code{\link{lmer}}, \code{\linkS4class{lmerenv}}}
\examples{
devcomp(lmer(Yield ~ 1 + (1|Batch), Dyestuff))
}
\keyword{models}
