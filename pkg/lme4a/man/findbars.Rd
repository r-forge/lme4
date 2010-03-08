\name{findbars}
\alias{findbars}
\alias{nobars}
\alias{subbars}
\title{Determine Random-Effects Expressions from a Formula}
\description{
  \code{findbars}:\cr
  From the right hand side of a formula for a mixed-effects model,
  determine the pairs of expressions that are separated by the
  vertical bar operator.

  \code{nobars}:\cr
  Omit terms separated by vertical bars in a formula, i.e.,
  remove the random-effects terms from a mixed-effects formula,
  thereby producing the fixed-effects formula.

  \code{subbars}:\cr
  \dQuote{Sub[stitute] Bars}:
  Substitute the \code{+} function for the \code{|} function in a
  mixed-model formula.  This provides a formula suitable for the current
  \code{\link{model.frame}} function.
}
\usage{
findbars(expr)
nobars(term)
subbars(term)
}
\arguments{
  \item{expr, term}{the right-hand side of a mixed-model formula}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
 a \code{\link{formula}}
}
%\author{Douglas Bates}
\seealso{
 \code{\link{formula}},
 \code{\link{model.frame}}, \code{\link{model.matrix}}.
}
\examples{
}
\keyword{models}
