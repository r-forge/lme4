\name{GaussQuad}
\alias{GaussQuad}
\title{
  Create knot and weight vectors for Gaussian quadrature
}
\description{
  Evaluate the positions (the \dQuote{knots}) and the weights for
  various Gaussian quadrature rules.  The quadrature rules are with
  respect to a kernel as described in the details section.  Optional
  arguments \code{a}, \code{b}, \code{alpha} and \code{beta} are used to
  generalize the rules.
}
\usage{
GaussQuad(n, rule, a = 0, b = 1, alpha = 0, beta = 0)
}
\arguments{
  \item{n}{integer - number of quadrature points}
  \item{rule}{character - a partial match to one of the quadrature rule
    names.  See the details section.}
  \item{a}{numeric scalar - an optional shift parameter or interval endpoint}
  \item{b}{numeric scalar - an optional scale parameter or interval endpoint}
  \item{alpha}{numeric scalar - an optional power}
  \item{beta}{numeric scalar - another optional power}
}
\details{
  The possible values for the \code{rule} character string and the
  corresponding integrals are:
  \describe{
    \item{\dQuote{Legendre}}{\eqn{\int_a^b f(x)\,dx}{int_a^b f(x) dx}}
    \item{\dQuote{Chebyshev}}{\eqn{\int_a^b f(x)((b-x)*(x-a))^{-0.5}\,dx}{int_a^b f(x) ((b-x)*(x-a))^(-0.5) dx}}
    \item{\dQuote{Gegenbauer}}{\eqn{\int_a^b f(x)((b-x)*(x-a))^\alpha\,dx}{int_a^b f(x) ((b-x)*(x-a))^alpha dx}}
    \item{\dQuote{Jacobi}}{\eqn{\int_a^b f(x)(b-x)^\alpha (x-a)^\beta\,dx}{int_a^b f(x) (b-x)^alpha*(x-a)^beta}}
    \item{\dQuote{Laguerre}}{\eqn{\int_a^\infty f(x)(x-a)^\alpha*exp(-b*(x-a))\,dx}{int_a^inf f(x)(x-a)^alpha*exp(-b*(x-a))}}
    \item{\dQuote{Hermite}}{\eqn{\int_{-\infty}^\infty f(x)|x-a|^\alpha \exp(-b*(x-a)^2)\,dx}{int_-inf^inf f(x)|x-a|^alpha*exp(-b*(x-a)^2) dx}}
    \item{\dQuote{Exponential}}{\eqn{\int_a^b f(x)|x-(a+b)/2.0|^\alpha\,dx}{int_-inf^inf f(x)|x-(a+b)/2.0|^alpha dx}}
    \item{\dQuote{Rational}}{\eqn{\int_a^\infty f(x)(x-a)^\alpha*(x+b)^\beta\,dx}{int_a^inf f(x)(x-a)^alpha*(x+b)^beta}}
  }
}
\value{
  A list with components
  \item{knots}{The positions at which to evaluate the function to be integrated}
  \item{weights}{The numeric weights to be applied.}
}
\references{
  The original FORTRAN implementation is from
  Sylvan Elhay, Jaroslav Kautsky (1987),
  \dQuote{Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature}, \bold{ACM Transactions on Mathematical Software},
  Volume \bold{13}, Number 4, December 1987, pages 399--415.

  One of the Eispack routines, an implicit QL algorithms from
  Roger Martin and James Wilkinson (1968), \dQuote{The Implicit QL Algorithm},
  \bold{Numerische Mathematik}, Volume \bold{12}, Number 5, December
  1968, pages 377--383, is used in a modified form.

  The C++ implementation by John Burkardt was modified to a C++ class by
  Douglas Bates.
}
%\author{}
%\note{}
%\seealso{}
\examples{
do.call(data.frame, GaussQuad(5L, "Hermite")) ## 5-point "physicist" Gauss-Hermite rule
do.call(data.frame, GaussQuad(5L, "H", b=0.5)) ## 5-point "probabilist" Gauss-Hermite rule
}
\keyword{math}
