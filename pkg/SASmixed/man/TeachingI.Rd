\name{TeachingI}
\alias{TeachingI}
\title{Teaching Methods I}
\description{
    The \code{TeachingI} data frame has 96 rows and 7 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{Method}{
	    a factor with levels \code{1} to \code{3}
	}
	\item{Teacher}{
	    a factor with levels \code{1} to \code{4}
	}
	\item{Gender}{
	    a factor with levels \code{f} and \code{m} 
	}
	\item{Student}{
	    a factor with levels \code{1} to \code{4}
	}
	\item{score}{
	    a numeric vector
	}
	\item{Experience}{
	    a numeric vector
	}
	\item{uTeacher}{
	    an ordered factor with levels
	}
    }
}
%\details{}
\source{
    Littel, R. C., Milliken, G. A., Stroup, W. W., and Wolfinger,
    R. D. (1996), \emph{SAS System for Mixed Models}, SAS Institute
    (Data Set 5.6).
}
\examples{
str(TeachingI)
}
\keyword{datasets}
