\name{Animal}
\alias{Animal}
\non_function{}
\title{Animal breeding experiment}
\description{
    The \code{Animal} data frame has 20 rows and 3 columns giving the
    average daily weight gains for animals with different genetic
    backgrounds.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{Sire}{
	    a factor denoting the sire. (5 levels)
	}
	\item{Dam}{
	    a factor denoting the dam. (2 levels)
	}
	\item{AvgDailyGain}{
	    a numeric vector of average daily weight gains
	}
    }
}\details{
    This appears to be a constructed data set.
}
\source{
    Littel, R. C., Milliken, G. A., Stroup, W. W., and Wolfinger,
    R. D. (1996), \emph{SAS System for Mixed Models}, SAS Institute
    (Data Set 6.4).
}
\examples{
str(Animal)
}
\keyword{datasets}
