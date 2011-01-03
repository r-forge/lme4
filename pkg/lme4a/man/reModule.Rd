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
    \item{\code{L}}{\code{\linkS4class{CHMfactor}} - sparse Cholesky
      factor (read only)} 
    \item{\code{Lambda}}{\code{\linkS4class{dgCMatrix}} - relative
      covariance factor (read only)} 
    \item{\code{Lambdap}}{external pointer to \code{Lambda} (read only)}
    \item{\code{Lind}}{integer - 1-based index vector into theta for
      \code{Lambda@x} (read only)}
    \item{\code{Lp}}{external pointer to \code{L} (read only)}
    \item{\code{Ut}}{external pointer to weighted, orthogonal design
      matrix (read only)}
    \item{\code{Zt}}{\code{\linkS4class{dgCMatrix}} - transpose of the
      model matrix for the random effects (read only)}
    \item{\code{b}}{numeric - random effects on original scale}
    \item{\code{cu}}{numeric - intermediate solution for u (read only)}
    \item{\code{incr}}{numeric - increment for \code{u}}
    \item{\code{ldL2}}{numeric - logarithm of the square of the
      determinant of \code{L} (read only)}
    \item{\code{linPred}}{numeric - linear predictor contribution based
      on current \code{u} (read only)}
    \item{\code{lower}}{numeric - lower bounds on components of
      \code{theta}}
    \item{\code{sqrLenU}}{numeric - squared length of u, the orthogonal
      random effects (read only)}
    \item{\code{theta}}{numeric - variance component parameters}
    \item{\code{u}}{numeric - orthogonal random effects vector (read
      only)}
    \item{\code{u0}}{numeric - base orthogonal random effects vector
      (read only)}
  }
}
\section{Methods}{
  \describe{
    \item{\code{linPred1(fac)}}{update \code{u} as \code{u0 + fac*incr}
      and return linear predictor contribution}
    \item{\code{reweight(sqrtXwt, wtres, useU0)}}{update L, Ut and cu
      for new weights.  The value of \code{u} is taken from \code{u0} or
      from \code{u0} according to whether or not \code{useU0} is TRUE}
    \item{\code{setU(u, delta, factor)}}{set a new value of u, possibly
      with an increment and step.}
    \item{\code{solveU()}}{solve for \code{u} (or the increment for
      \code{u}) only. Returns squared length of c1}
    \item{\code{updateIncr(cu)}}{solve for \code{incr} given the updated
      \code{cu} from the feModule's \code{updateBeta} method}
    \item{\code{updateU(cu)}}{solve for \code{u} given the updated
      \code{cu} from the feModule's \code{updateBeta} method}
  }
}
%\details{}
%\source{}
%\references{}
%\examples{}
\keyword{datasets}
