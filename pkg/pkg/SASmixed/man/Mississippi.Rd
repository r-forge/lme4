\name{Mississippi}
\alias{Mississippi}
\non_function{}
\title{Nitrogen concentrations in the Mississippi River}
\description{
The \code{Mississippi} data frame has 37 rows and 3 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{influent}{
	    an ordered factor with levels
	    \code{3} < \code{5} < \code{2} < \code{1} < \code{4} < \code{6}
	}
	\item{y}{
	    a numeric vector
	}
	\item{Type}{
	    a factor with levels
	    \code{1} 
	    \code{2} 
	    \code{3} 
	}
    }
}
\source{
    Littel, R. C., Milliken, G. A., Stroup, W. W., and Wolfinger,
    R. D. (1996), \emph{SAS System for Mixed Models}, SAS Institute
    (Data Set 4.2).
}
\examples{
options(contrasts = c(unordered = "contr.SAS", ordered = "contr.poly"))
str(Mississippi)
#plot(Mississippi)
fm1Miss <- lmer(y ~ 1 + (1|influent), Mississippi)
summary(fm1Miss)        # compare with output 4.1, p. 142
fm1MLMiss <- update(fm1Miss, method = "ML")
summary(fm1MLMiss)        # compare with output 4.2, p. 143
ranef(fm1MLMiss)          # BLUP's of random effects on p. 144
ranef(fm1Miss)            # BLUP's of random effects on p. 142
#intervals(fm1Miss)       # interval estimates of variance components
fm2Miss <- lmer(y ~ Type+(1|influent), Mississippi, method = "REML")
summary(fm2Miss)         # compare to output 4.8 and 4.9, pp. 150-152
anova(fm2Miss)
}
\keyword{datasets}
