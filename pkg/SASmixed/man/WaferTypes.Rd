\name{WaferTypes}
\alias{WaferTypes}
\non_function{}
\title{Data on different types of silicon wafers}
\description{
    The \code{WaferTypes} data frame has 144 rows and 8 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{Group}{
	    a factor with levels \code{1} to \code{4}
	}
	\item{Temperature}{
	    an ordered factor with levels
	    \code{900} < \code{1000} < \code{1100}
	}
	\item{Type}{
	    a factor with levels \code{A} and \code{B} 
	}
	\item{Wafer}{
	    a numeric vector
	}
	\item{Site}{
	    a numeric vector
	}
	\item{delta}{
	    a numeric vector
	}
	\item{Thick}{
	    a numeric vector
	}
	\item{uWafer}{
	    an ordered factor giving a unique code to each group,
	    temperature, type and wafer combination.
	}
    }
}
%\details{}
\source{
    Littel, R. C., Milliken, G. A., Stroup, W. W., and Wolfinger,
    R. D. (1996), \emph{SAS System for Mixed Models}, SAS Institute
    (Data Set 5.8).
}
\examples{
str(WaferTypes)
}
\keyword{datasets}
