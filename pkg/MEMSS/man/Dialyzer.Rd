\name{Dialyzer}
\alias{Dialyzer}
\title{High-Flux Hemodialyzer}
\description{
    The \code{Dialyzer} data frame has 140 rows and 5 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{Subject}{
	    a factor with levels \code{A} to \code{T}
	}
	\item{QB}{
	    a factor with levels
	    \code{200} and 
	    \code{300} giving the bovine blood flow rate (dL/min).
	}
	\item{pressure}{
	    the transmembrane pressure (dmHg).
	}
	\item{rate}{
	    the hemodialyzer ultrafiltration rate (mL/hr).
	}
	\item{index}{
	    index of observation within subject---1 through 7.
	}
    }
}
\details{
    Vonesh and Carter (1992) describe data measured on high-flux
    hemodialyzers to assess their \emph{in vivo} ultrafiltration
    characteristics. The ultrafiltration rates (in mL/hr) of 20 high-flux
    dialyzers were measured at seven different transmembrane pressures (in
    dmHg). The \emph{in vitro} evaluation of the dialyzers used bovine
    blood at flow rates of either 200~dl/min or 300~dl/min. The data,
    are also analyzed in Littell, Milliken, Stroup, and Wolfinger (1996).
}
\source{
    Pinheiro, J. C. and Bates, D. M. (2000), \emph{Mixed-Effects Models in S
	and S-PLUS}, Springer, New York. (Appendix A.6)

    Vonesh, E. F. and Carter, R. L. (1992), Mixed-effects nonlinear
    regression for unbalanced repeated measures, \emph{Biometrics},
    \bold{48}, 1-18.

    Littell, R. C., Milliken, G. A., Stroup, W. W. and Wolfinger,
    R. D. (1996), \emph{SAS System for Mixed Models}, SAS Institute,
    Cary, NC.
}
\examples{
str(Dialyzer)
}
\keyword{datasets}
