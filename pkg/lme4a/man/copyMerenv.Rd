\name{copyMerenv}
\alias{copyMerenv}
\title{Create a Deep Copy of a 'merenv' Object}
\description{
  As objects of class \code{\linkS4class{merenv}}, mainly contain an
  environment, the standard \dQuote{copying} \code{y <- x} call does not
  effectively create a copy of \code{x}, but rather a \dQuote{shallow}
  copy with the peculiar property that changing the objects in
  \code{\link{env}(x)}, those in \code{\link{env}(y)} are automatically
  changed along, as they \emph{are} effectively the same objects.

  Instead, \code{y <- copyMerenv(x)} makes \code{y} a \dQuote{deep} copy
  of \code{x}.
}
\usage{
copyMerenv(x)
}
\arguments{
  \item{x}{an object of class \code{\linkS4class{merenv}}.}
}
\value{
  an object \code{\link{identical}} to \code{x}, but with objects as
  content that are true \emph{copies} of those in \code{x}.
}
\seealso{
  \code{\link{lmer}}, \code{\linkS4class{merenv}}.
}
%\examples{
%}
\keyword{utilities}
