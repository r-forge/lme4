\name{Cultivation}
\alias{Cultivation}
\non_function{}
\title{Bacterial innoculation applied to grass cultivars}
\description{
    The \code{Cultivation} data frame has 24 rows and 4 columns of data
    from an experiment on the effect on dry weight yield of three
    bacterial inoculation treatments applied to two grass cultivars.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{Block}{
	    a factor with levels \code{1} to \code{4}
	}
	\item{Cult}{
	    the cultivar factor with levels \code{a} and \code{b}
	}
	\item{Inoc}{
	    the innoculant factor with levels \code{con}, \code{dea} and
	    \code{liv}
	}
	\item{drywt}{
	    a numeric vector of dry weight yields
	}
    }
}
\source{
    Littel, R. C., Milliken, G. A., Stroup, W. W., and Wolfinger,
    R. D. (1996), \emph{SAS System for Mixed Models}, SAS Institute
    (Data Set 2.2(a)).

    Littel, R. C., Freund, R. J., and Spector, P. C. (1991), \emph{SAS
	System for Linear Models, Third Ed.}, SAS Institute.
}
\examples{
options(contrasts = c(unordered = "contr.SAS", ordered = "contr.poly"))
str(Cultivation)
xtabs(~Block+Cult, Cultivation)
fm1Cult <- lmer(drywt ~ Inoc * Cult + (1|Block) + (1|Cult), Cultivation)
summary(fm1Cult)   # compare with Output 2.10, page 58
anova(fm1Cult)
fm2Cult <- lmer(drywt ~ Inoc + Cult + (1|Block) + (1|Cult), Cultivation)
anova(fm2Cult)
fm3Cult <- lmer(drywt ~ Inoc + (1|Block) + (1|Cult), Cultivation)
anova(fm3Cult)
summary(fm3Cult)
}
\keyword{datasets}
