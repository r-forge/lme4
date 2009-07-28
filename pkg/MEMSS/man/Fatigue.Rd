\name{Fatigue}
\alias{Fatigue}
\title{Cracks caused by metal fatigue}
\description{
  The \code{Fatigue} data frame has 262 rows and 3 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Path}{
      the test path (or test unit) identifier - a factor with levels
      \code{A} to \code{U}.
    }
    \item{cycles}{
      number of test cycles at which the measurement is made
      (millions of cycles).
    }
    \item{relLength}{
      relative crack length (dimensionless).
    }
  }
}
\details{
  These data are given in Lu and Meeker (1993) where they state
  ``We obtained the data in Table 1 visually from figure 4.5.2 on page 242
  of Bogdanoff and Kozin (1985).''
  The data represent the growth of cracks in metal for 21 test units.
  An initial notch of length 0.90 inches was made on each unit which then
  was subjected to several thousand test cycles.  After every 10,000
  test cycles the crack length was measured.  Testing was stopped if the
  crack length exceeded 1.60 inches, defined as a failure, or at 120,000
  cycles.
}
\source{
  Lu, C. Joseph , and Meeker, William Q. (1993), Using degradation
  measures to estimate a time-to-failure distribution,
  \emph{Technometrics}, \bold{35}, 161-174 
}
\examples{
str(Fatigue)
xyplot(relLength ~ cycles | Path, Fatigue, type = c("g", "b"),
       aspect = 'xy', xlab = "Number of test cycles (millions)",
       ylab = "Relative crack length (dimensionless)",
       layout = c(7,3))
}
\keyword{datasets}
