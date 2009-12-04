\name{sigma}
\title{Estimate of common scale parameter}
\usage{
sigma(object, ...)
}
\alias{sigma}
\alias{sigma,mer-method}
\alias{sigma,lmerenv-method}
\description{
  Return the estimate of the common scale parameter from a fitted model
  object, typically a mixed-effects model of class
  \code{\linkS4class{lmerenv}}.
}
\arguments{
  \item{object}{a fitted model object of a suitable class}
  \item{\dots}{optional additional arguments.  None are used at present.}
}
\details{
  Methods for this function return the estimate of the standard
  deviation of the residual or \dQuote{per-observation} noise term in a
  mixed-effects model.  For linear mixed effects models the REML
  estimate is returned if the model was fit according to the REML
  criterion, otherwise the maximum likelihood estimate is returned.

  For generalized linear mixed models that do not incorporate a common
  scale parameter, in particular those for which \code{family =
  binomial} or \code{family = poisson}, \code{NA} is returned.
}
\value{A non-negative numeric scalar}
\seealso{\code{\link{lmer}}, \code{\link{glmer}}}
\examples{
sigma(fm1 <- lmer(Yield ~ 1 + (1|Batch), Dyestuff)) # REML estimate
sigma(update(fm1, REML = FALSE))     # ML estimate
}
\keyword{models}
\keyword{methods}
