\name{Exam}
\alias{Exam}
\docType{data}
\title{Exam scores from inner London}
\description{
  Exam scores of 4,059 students from 65 schools in Inner London.
}
\usage{data(Exam)}
\format{
  A data frame with 4059 observations on the following 9 variables.
  \describe{
    \item{school}{School ID - a factor.}
    \item{normexam}{Normalized exam score.}
    \item{schgend}{School gender - a factor.  Levels are \code{mixed},
      \code{boys}, and \code{girls}.}
    \item{schavg}{School average of intake score.}
    \item{vr}{Student level Verbal Reasoning (VR) score band at intake -
      a factor.  Levels are \code{bottom 25\%}, \code{mid 50\%}, and
      \code{top 25\%}.}
    \item{intake}{Band of student's intake score - a factor.
      Levels are \code{bottom 25\%}, \code{mid 50\%} and \code{top
	25\%}./}
    \item{standLRT}{Standardised LR test score.}
    \item{sex}{Sex of the student - levels are \code{F} and \code{M}.}
    \item{type}{School type - levels are \code{Mxd} and \code{Sngl}.}
    \item{student}{Student id (within school) - a factor}
  }
}
%\details{}
\source{
  \url{http://multilevel.ioe.ac.uk/softrev/exam.html}
}
\references{
  Goldstein, H., Rasbash, J., et al (1993). A multilevel analysis of
  school examination results. \emph{Oxford Review of Education} 19: 425-433 
}
\examples{
str(Exam)
summary(Exam)
(fm1 <- lmer(normexam ~ standLRT + sex + schgend + (1|school), Exam))
(fm2 <- lmer(normexam ~ standLRT*sex + schgend + (1|school), Exam))
(fm3 <- lmer(normexam ~ standLRT*sex + schgend + (1|school), Exam))
}
\keyword{datasets}
