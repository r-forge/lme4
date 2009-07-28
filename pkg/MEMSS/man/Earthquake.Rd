\name{Earthquake}
\alias{Earthquake}
\title{Earthquake Intensity}
\description{
    The \code{Earthquake} data frame has 182 rows and 5 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{Quake}{
	    a factor with levels \code{A} to \code{U}
	}
	\item{Richter}{
	    the intensity of the earthquake on the Richter scale
	}
	\item{distance}{
	    the distance from the seismological measuring station to the
	    epicenter of the earthquake (km)
	}
	\item{soil}{
	    a factor with levels \code{S} (soil) and \code{R} (rock)
	    giving the soil condition at the measuring station 
	}
	\item{accel}{
	    maximum horizontal acceleration observed (g).
	}
    }
}
\details{
    Measurements recorded at available seismometer locations for 23
    large earthquakes in western North America between 1940 and 1980.
    They were originally given in Joyner and Boore (1981); are
    mentioned in Brillinger (1987); and are analyzed in Davidian
    and Giltinan (1995).
}
\source{
    Pinheiro, J. C. and Bates, D. M. (2000), \emph{Mixed-Effects Models in S
	and S-PLUS}, Springer, New York.  (Appendix A.8)

    Davidian, M. and Giltinan, D. M. (1995), \emph{Nonlinear Models for
	Repeated Measurement Data},  Chapman and Hall, London.

    Joyner and Boor (1981), Peak horizontal acceleration and velocity from
    strong-motion records including records from the 1979 Imperial
    Valley, California, earthquake, \emph{Bulletin of the Seismological
	Society of America}, \bold{71}, 2011-2038.

    Brillinger, D. (1987), Comment on a paper by C. R. Rao,
    \emph{Statistical Science}, \bold{2}, 448-450.
}
\examples{
str(Earthquake)
}
\keyword{datasets}
