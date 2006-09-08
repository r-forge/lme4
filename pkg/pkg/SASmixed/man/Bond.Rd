\name{Bond}
\alias{Bond}
\non_function{}
\title{Strengths of metal bonds}
\description{
    The \code{Bond} data frame has 21 rows and 3 columns of data on the
    strength required to break metal bonds according to the metal and
    the ingot.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{pressure}{
	    a numeric vector of pressures required to break the bond
	}
	\item{Metal}{
	    a factor with levels \code{c}, \code{i} and \code{n}
	    indicating the metal involved (copper, iron or nickel).
	}
	\item{Ingot}{
	    an ordered factor indicating the ingot of the composition
	    material.
	}
    }
}
\source{
    Littel, R. C., Milliken, G. A., Stroup, W. W., and Wolfinger,
    R. D. (1996), \emph{SAS System for Mixed Models}, SAS Institute
    (Data Set 1.2.4).

    Mendenhall, M., Wackerly, D. D. and Schaeffer, R. L. (1990),
    \emph{Mathematical Statistics}, Wadsworth (Exercise 13.36).
}
\examples{
options(contrasts = c(unordered = "contr.SAS", ordered = "contr.poly"))
fm1Bond <- lmer(pressure ~ Metal + (1|Ingot), Bond)
summary(fm1Bond)               # compare with output 1.1 on p. 6
anova(fm1Bond)
}
\keyword{datasets}
