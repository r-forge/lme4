\name{BIB}
\alias{BIB}
\non_function{}
\title{Data from a balanced incomplete block design}
\description{
    The \code{BIB} data frame has 24 rows and 5 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{Block}{
	    an ordered factor with levels
	    \code{1} < \code{2} < \code{3} < \code{8} < \code{5} < \code{4} < \code{6} < \code{7}
	}
	\item{Treatment}{
	    a treatment factor with levels \code{1} to \code{4}.
	}
	\item{y}{
	    a numeric vector representing the response
	}
	\item{x}{
	    a numeric vector representing the covariate
	}
	\item{Grp}{
	    a factor with levels
	    \code{13} and 
	    \code{24} 
	}
    }
}
\details{
    These appear to be constructed data.
}
\source{
    Littel, R. C., Milliken, G. A., Stroup, W. W., and Wolfinger,
    R. D. (1996), \emph{SAS System for Mixed Models}, SAS Institute
    (Data Set 5.4).
}
\examples{
options(contrasts = c(unordered = "contr.SAS", ordered = "contr.poly"))
if (require("lattice", quietly = TRUE, character = TRUE)) {
  xyplot(y ~ x | Block, BIB, groups = Treatment, type = c("g", "p"),
         aspect = "xy", auto.key = list(points = TRUE, space = "right",
         lines = FALSE))
}
fm1BIB <- lmer(y ~ Treatment * x + (1 | Block), BIB)
summary(fm1BIB)   # compare with Output 5.7, p. 188
anova(fm1BIB)     # strong evidence of different slopes
fm2BIB <- lmer(y ~ Treatment + x : Grp + (1 | Block), BIB)
summary(fm2BIB)   # compare with Output 5.9, p. 193
anova(fm2BIB)
}
\keyword{datasets}
