% $Id: pooledSD.Rd,v 1.2 2003/11/01 04:13:19 deepayan Exp $
\name{pooledSD}
\title{Extract pooled standard deviation}
\usage{
pooledSD(x, \dots)
}
\alias{pooledSD}
\arguments{
  \item{x}{an object inheriting from class \code{lmList}.}
  \item{\dots}{optional additional arguments (none used at present)}  
}
\description{
  The pooled estimated standard deviation is obtained by adding together
  the residual sum of squares for each non-null element of
  \code{x}, dividing by the sum of the corresponding residual
  degrees-of-freedom, and taking the square-root.
}
\value{
  the pooled standard deviation for the non-null elements of
  \code{x}, with an attribute \code{df} with the number of
  degrees-of-freedom used in the estimation.
}
\seealso{\code{\link{lmList}}, \code{\link{lm}}}
\examples{
pooledSD(lmList(weight ~ Time | Chick, do.call("data.frame", ChickWeight)))
}
\keyword{models}
