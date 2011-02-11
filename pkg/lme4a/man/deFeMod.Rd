\name{deFeMod}
\alias{deFeMod}
\alias{coerce,Rcpp_deFeMod,deFeMod-method}
\alias{show,Rcpp_deFeMod-method}
\docType{class}
\title{
  An C++ class for information related to the fixed-effects parameters.
}
\description{
  This C++Class encapsulates the information associated with the fixed
  effects in a mixed-effects model.
}
\section{Constructors}{
  \describe{
    \item{\code{deFeMod(Rcpp::S4, int)}:}{from an S4
      \code{\linkS4class{deFeMod}} object and the number of
      observations, \code{n}}
    \item{\code{deFeMod(X, n, q)}}{from the dense model matrix, \code{X}
      and dimensions}
  }
}
\section{Fields}{
  \describe{
    \item{\code{RX}}{\code{\linkS4class{Cholesky}} - dense Cholesky
      factor (read only)}
    \item{\code{RZX}}{\code{\linkS4class{dgeMatrix}} - Off-diagonal
      block in large Cholesky factor (read only)}
    \item{\code{UtV}}{\code{\linkS4class{dgeMatrix}} - Weighted
      crossproduct (read only)}
    \item{\code{V}}{\code{\linkS4class{dgeMatrix}} - Scaled model matrix
       (read only)}
    \item{\code{VtV}}{\code{\linkS4class{dpoMatrix}} - Weighted
      crossproduct (read only)}
    \item{\code{Vtr}}{numeric - Weighted crossproduct of model matrix
      and residuals (read only)}
    \item{\code{X}}{\code{\linkS4class{dgeMatrix}} - original model
      matrix (read only)}
    \item{\code{coef}}{numeric - coefficient vector (read only)}
    \item{\code{coef0}}{numeric - base coefficient vector (read only)}
    \item{\code{incr}}{numeric - increment for coefficient vector (read
      only)}
    \item{\code{ldRX2}}{numeric - log of the square of the determinant
      of RX (read only)}
    \item{\code{linPred}}{numeric - linear predictor contribution (read
      only)}
  }
}
\section{Methods}{
  \describe{
    \item{\code{linPred1(fac)}}{numeric - update \code{coef} to
      \code{coef0 + fac*incr} and return linear predictor contribution}
    \item{\code{reweight(sqrtXwt, wtres)}}{update V and Vtr for new
      weights.}
    \item{\code{updateBeta(cu)}}{update \code{coef} given 
      \code{cu}}
    \item{\code{updateIncr(cu)}}{update \code{incr} given 
      \code{cu}}
    \item{\code{updateRzxRx(Lambda, L)}}{update update the triangular
      factor sections given \code{Lambda} and \code{L}}
    \item{\code{updateRzxpRxp(Lambdap, Lp)}}{update the
      triangular factor sections given external pointers to
      \code{Lambda} and \code{L}}
    \item{\code{updateUtV(Utp)}}{update \code{UtV} 
      given an external pointers to
      \code{Ut}}
  }
}
%\details{}
%\source{}
%\references{}
%\examples{}
\keyword{datasets}
