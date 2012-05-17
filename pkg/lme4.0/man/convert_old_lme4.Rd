\name{convert_old_lme4}
\alias{convert_old_lme4}
\title{
  Convert old lme4 objects to lme4.0
}
\description{
  Given a fitted model created by the old (<= 0.999375-42) \code{lme4} package,
  convert it to be usable with the \code{lme4.0} package
}
\usage{
convert_old_lme4(x)
}
\arguments{
  \item{x}{
    A fitted model from \code{lme4} (of class \code{mer})
}
}
\details{
  The function simply manipulates the \code{package} attribute of
  the object's class to allow \code{lme4.0} to recognize it.
}
\value{
  A \code{mer} object usable by \code{lme4.0}
}
\author{
  Ben Bolker
}
\examples{
## an example file containing old mer objects
(L <- load(system.file("old_lme4_save.RData",package="lme4.0")))
convert_old_lme4(fm1)
convert_old_lme4(gm1)
}
\keyword{misc}

