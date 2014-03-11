\name{Glucose}
\alias{Glucose}
\title{Glucose levels over time}
\description{
    The \code{Glucose} data frame has 378 rows and 4 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{Subject}{
	    a factor with levels \code{A} to \code{F}
	}
	\item{Time}{
	    a numeric vector
	}
	\item{conc}{
	    a numeric vector of glucose levels
	}
	\item{Meal}{
	    an ordered factor with levels
	    \code{2am} < \code{6am} < \code{10am} < \code{2pm} < \code{6pm} < \code{10pm}
	}
    }
}
\source{
    Hand, D. and Crowder, M. (1996), \emph{Practical Longitudinal Data
	Analysis}, Chapman and Hall, London.
}
\examples{
require(lattice)
str(Glucose)
xyplot(conc ~ Time | Meal * Subject, Glucose)
}
\keyword{datasets}
