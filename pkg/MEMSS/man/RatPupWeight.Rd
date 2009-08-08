\name{RatPupWeight}
\alias{RatPupWeight}
\title{The weight of rat pups}
\description{
    The \code{RatPupWeight} data frame has 322 rows and 5 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{weight}{
	    a numeric vector
	}
	\item{sex}{
	    a factor with levels
	    \code{Male} 
	    \code{Female} 
	}
	\item{Litter}{an ordered factor}
	\item{Lsize}{a numeric vector}
	\item{Treatment}{
	    an ordered factor with levels
	    \code{Control} < \code{Low} < \code{High}
	}
    }
}
\source{
    Pinheiro, J. C. and Bates, D. M. (2000), \emph{Mixed-Effects Models in S
	and S-PLUS}, Springer, New York.
}
\examples{
str(RatPupWeight)
}
\keyword{datasets}
