\name{VerbAgg}
\alias{VerbAgg}
\docType{data}
\title{Verbal Aggression item responses}
\description{
  These are the item responses to a questionaire on verbal aggression.
  These data are used throughout De Boeck and Wilson,
  \emph{Explanatory Item Response Models} (Springer, 2004) to illustrate
  various forms of item response models.
}
%\usage{data(VerbAgg)}
\format{
  A data frame with 7584 observations on the following 13 variables.
  \describe{
    \item{\code{Anger}}{the subject's Trait Anger score as measured on
      the State-Trait Anger Expression Inventory (STAXI)}
    \item{\code{Gender}}{the subject's gender - a factor with levels
      \code{M} and \code{F}}
    \item{\code{item}}{the item on the questionaire, as a factor}
    \item{\code{resp}}{the subject's response to the item - an ordered
      factor with levels \code{no} < \code{perhaps} < \code{yes}}
    \item{\code{id}}{the subject identifier, as a factor}
    \item{\code{btype}}{behavior type - a factor with levels \code{curse} \code{scold} \code{shout}}
    \item{\code{situ}}{situation type - a factor with levels \code{other} \code{self}}
    \item{\code{mode}}{behavior mode - a factor with levels \code{want} \code{do}}
    \item{\code{r2}}{dichotomous version of the response - a factor with levels \code{N} \code{Y}}
    \item{\code{Express}}{a numeric vector of one of the contrasts for
      \code{btype} used in De Boeck and Wilson (2004)}
    \item{\code{Blame}}{a numeric vector of the other contrast for
      \code{btype} used in De Boeck and Wilson (2004)}
    \item{\code{OtherSelf}}{a numeric vector of a contrast for the \code{situ}}
    \item{\code{DoWant}}{a numeric vector of a contrast for the \code{mode}}
  }
}
%\details{}
\source{
  \url{http://bear.soe.berkeley.edu/EIRM/}
}
\references{
  De Boeck and Wilson (2004), \emph{Explanatory Item Response Models}, Springer.
}
\examples{
str(VerbAgg)
xtabs(~ item + resp, VerbAgg)
xtabs(~ btype + resp, VerbAgg)
round(100 * ftable(prop.table(xtabs(~ situ + mode + resp, VerbAgg), 1:2), 1))
person <- unique(subset(VerbAgg, select = c(id, Gender, Anger)))
if (require(lattice)) {
    densityplot(~ Anger, person, groups = Gender, auto.key = list(columns = 2),
                xlab = "Trait Anger score (STAXI)")
}
}
\keyword{datasets}
