\name{HR}
\alias{HR}
\non_function{}
\title{Heart rates of patients on different drug treatments}
\description{
    The \code{HR} data frame has 120 rows and 5 columns of the heart
    rates of patients under one of three possible drug treatments.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{Patient}{
	    an ordered factor indicating the patient.
	}
	\item{Drug}{
	    the drug treatment - a factor with levels \code{a},
	    \code{b} and \code{p} where \code{p} represents the placebo.
	}
	\item{baseHR}{
	    the patient's base heart rate
	}
	\item{HR}{
	    the observed heart rate at different times in the experiment
	}
	\item{Time}{
	    the time of the observation
	}
    }
}
\source{
    Littel, R. C., Milliken, G. A., Stroup, W. W., and Wolfinger,
    R. D. (1996), \emph{SAS System for Mixed Models}, SAS Institute
    (Data Set 3.5).
}
\examples{
options(contrasts = c(unordered = "contr.SAS", ordered = "contr.poly"))
if (require("lattice", quietly = TRUE, character = TRUE)) {
  xyplot(HR ~ Time | Patient, HR, type = c("g", "p", "r"), aspect = "xy",
         index.cond = function(x, y) coef(lm(y ~ x))[1],
         ylab = "Heart rate (beats/min)")
}
(fm1HR <- lmer(HR ~ Time * Drug + baseHR + (Time|Patient), HR))  # linear trend in time
anova(fm1HR)
\dontrun{
fm2HR <- update(fm1HR, weights = varPower(0.5)) # use power-of-mean variance
summary(fm2HR)
intervals(fm2HR)             # variance function does not seem significant
anova(fm1HR, fm2HR)         # confirm with likelihood ratio
}
(fm3HR <- lmer(HR ~ Time + Drug + baseHR + (Time|Patient), HR))
anova(fm3HR)
(fm4HR <- lmer(HR ~ Time + baseHR + (Time|Patient), HR))  # remove Drug term
anova(fm4HR)
}
\keyword{datasets}
