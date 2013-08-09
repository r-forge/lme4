\name{convert_old_lme4}
\alias{convert_old_lme4}
\title{Convert Old-lme4 Objects to lme4.0 Ones}
\description{
  Given a fitted model created by the old (<= 0.999999-2) \pkg{lme4} package,
  convert it to be usable with the \pkg{lme4.0} package.
}
\usage{
convert_old_lme4(x)
}
\arguments{
  \item{x}{fitted model from \pkg{lme4} (of class \code{\linkS4class{mer}}).}
}
\details{
  The function simply manipulates the \code{package} attribute of
  the object's \code{\link{class}} to allow \pkg{lme4.0} to recognize it.
}
\value{
  A \code{\linkS4class{mer}} object usable by \pkg{lme4.0}.
}
\author{
  Ben Bolker
}
\examples{
## an example file containing old mer objects
(load(system.file("old_lme4_save.RData", package="lme4.0")))
validObject(fm2 <- convert_old_lme4(fm1)); fm2
validObject(gm2 <- convert_old_lme4(gm1)); gm2
## they differ by package attribute only :
stopifnot(all.equal(fm1, fm2, check.attributes=FALSE, tol=0))
}
\keyword{misc}

