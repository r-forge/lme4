\name{Genetics}
\alias{Genetics}
\non_function{}
\title{Heritability data}
\description{
The \code{Genetics} data frame has 60 rows and 4 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{Location}{
	    a factor with levels \code{1} to \code{4} 
	}
	\item{Block}{
	    a factor with levels \code{1} to \code{3}
	}
	\item{Family}{
	    a factor with levels \code{1} to \code{5}
	}
	\item{Yield}{
	    a numeric vector of crop yields
	}
    }
}
\source{
    Littel, R. C., Milliken, G. A., Stroup, W. W., and Wolfinger,
    R. D. (1996), \emph{SAS System for Mixed Models}, SAS Institute
    (Data Set 4.5).
}
\examples{
str(Genetics)
if (require("lme4", quietly = TRUE, character = TRUE)) {
  options(contrasts = c(unordered = "contr.SAS", ordered = "contr.poly"))
  summary(fm1Gen <- lmer(Yield ~ Family + (1|Location/Block), Genetics))
}
}
\keyword{datasets}
