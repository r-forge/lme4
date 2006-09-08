\name{TeachingII}
\alias{TeachingII}
\non_function{}
\title{Teaching Methods II}
\description{
    The \code{TeachingII} data frame has 96 rows and 6 columns.
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
	\item{IQ}{
	    a numeric vector
	}
	\item{score}{
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
    (Data Set 5.7).
}
\examples{
str(TeachingII)
}
\keyword{datasets}
