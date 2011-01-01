\name{reModule}
\alias{reModule}
\docType{class}
\title{
  An C++ class for information related to random effects
}
\description{
  This C++Class encapsulates the information associated with the random
  effects in a mixed-effects model.
}
\section{Constructors}{
  \describe{
    \item{\code{reModule(Rcpp::S4)}:}{from an S4
      \code{\linkS4class{reModule}} object}
    \item{\code{reModule(Zt, Lambda, L, Lind, lower)}}{from components}
  }
}
\section{Fields}{
  \describe{
    \item{\code{L}}{chmFr - sparse Cholesky factor (read only)}
    \item{\code{Lambda}}{chmSp - relative covariance factor (read only)}
    \item{\code{Lambdap}}{external pointer to \code{Lambda} (read only)}
    \item{\code{Lind}}{integer - 1-based index vector into theta for
      \code{Lambda@x} (read only)}
    \item{\code{Lp}}{external pointer to \code{L} (read only)}
    \item{\code{Ut}}{external pointer to weighted, orthogonal design
      matrix (read only)}
    \item{\code{Zt}}{chmSp - transpose of the model matrix for the
      random effects (read only)}
    \item{\code{cu}}{numeric - intermediate solution for u (read only)}
    \item{\code{ldL2}}{numeric - logarithm of the square of the
      determinant of L (read only)}
    \item{\code{linPred}}{numeric - linear predictor contribution (read
      only)}
    \item{\code{sqrLenU}}{numeric - squared length of u, the orthogonal
      random effects (read only)}
    \item{\code{theta}}{numeric - variance component parameters}
    \item{\code{u}}{numeric - orthogonal random effects vector (read
      only)}
  }
}
\section{Methods}{
  \describe{
    \item{\code{reweight(sqrtXwt, wtres)}}{update L, Ut and cu for new
      weights.}
    \item{\code{setU(u, delta, factor)}}{set a new value of u, possibly
      with an increment and step.}
    \item{\code{solveU()}}{solve for \code{u} (or the increment for
      \code{u}) only. Returns squared length of c1}
    \item{\code{updateU(cu)}}{solve for \code{u} given the updated
      \code{cu} from the feModule's \code{updateBeta} method}
  }
}
%\details{}
%\source{}
%\references{}
%\examples{}
\keyword{datasets}
