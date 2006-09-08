\name{Semi2}
\alias{Semi2}
\non_function{}
\title{Oxide layer thicknesses on semiconductors}
\description{
    The \code{Semi2} data frame has 72 rows and 5 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{Source}{
	    a factor with levels \code{1} and \code{2} 
	}
	\item{Lot}{
	    a factor with levels \code{1} to \code{8}
	}
	\item{Wafer}{
	    a factor with levels \code{1} to \code{3}
	}
	\item{Site}{
	    a factor with levels \code{1} to \code{3}
	}
	\item{Thickness}{
	    a numeric vector
	}
    }
}
%\details{}
\source{
    Littel, R. C., Milliken, G. A., Stroup, W. W., and Wolfinger,
    R. D. (1996), \emph{SAS System for Mixed Models}, SAS Institute
    (Data Set 4.4).
}
\examples{
options(contrasts = c(unordered = "contr.SAS", ordered = "contr.poly"))
str(Semi2)
xtabs(~Lot + Wafer, Semi2)
Semi2$LotWafer <- with(Semi2, Lot:Wafer)
fm1Semi2 <- lmer(Thickness ~ 1 + (1|LotWafer) + (1|Lot), Semi2)
summary(fm1Semi2)       # compare with output 4.13, p. 156
fm2Semi2 <- lmer(Thickness ~ Source + (1|LotWafer) + (1|Lot), Semi2)
summary(fm2Semi2)       # compare with output 4.15, p. 159
anova(fm2Semi2)
Semi2$LotSrc <- with(Semi2, Lot:Source)
fm3Semi2 <- lmer(Thickness ~ Source + (1|LotWafer) + (1|LotSrc), Semi2)
summary(fm3Semi2)       # compare with output 4.17, p. 163
## This is not the same as the SAS model.
}
\keyword{datasets}
