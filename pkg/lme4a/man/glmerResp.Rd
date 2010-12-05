\name{glmerResp}
\alias{glmerResp}
\docType{class}
\title{
  An C++ class for \code{glmRespModule} objects
}
\description{
  The \code{"glmRespModule"} S4 class represents the information
  associated with the response in a \code{\link{glm}} or
  \code{\link{glmer}} model.  This C++Class encapsulates the
  information and provides methods accessible in \R{} or in \code{C++}. 
}
\section{Constructors}{
\describe{
  \item{\code{glmerResp(Rcpp::S4)}:}{from a
    \code{\linkS4class{glmRespMod}} object}
  \item{\code{glmerResp(Rcpp::List, Rcpp::NumericVector)}}{from a glm
    \code{\link{family}} and the response \code{y}}
  \item{\code{glmerResp(Rcpp::List, Rcpp::NumericVector, Rcpp::NumericVector)}}{from a glm
    \code{\link{family}}, the response \code{y} and the \code{weights}}
  \item{\code{glmerResp(Rcpp::List, Rcpp::NumericVector,
      Rcpp::NumericVector, Rcpp::NumericVector)}}{from a glm
    \code{\link{family}}, the response \code{y}, the \code{weights} and
    the \code{offset}}
  \item{\code{glmerResp(Rcpp::List, Rcpp::NumericVector,
      Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector)}}{from a glm
    \code{\link{family}}, the response \code{y}, the \code{weights},
    the \code{offset} and the number of replicates, \code{n} (only
    applies to the binomial family).}
}
}
\section{Fields}{
  \describe{
    \item{\code{devResid}}{numeric - deviance residuals (read only)}
    \item{\code{eta}}{numeric - linear predictor (read only)}
    \item{\code{family}}{character - name of glm family (read only)}
    \item{\code{link}}{character - name of link function (read only)}
    \item{\code{mu}}{numeric - mean vector (read only)}
    \item{\code{muEta}}{numeric - derivative of \code{mu} with respect
      to \code{eta} (readonly)}
    \item{\code{offset}}{numeric - the offset vector, always present
      even if zero (read only)}    
    \item{\code{residDeviance}}{numeric scalar - sum of squared deviance
      residuals (read only)}
    \item{\code{sqrtWrkWt}}{numeric - square root of the weights applied
      to the working residuals or working response (read only)}
    \item{\code{sqrtXwt}}{numeric - square root of the weights applied
      to the model matrix (read only)}
    \item{\code{sqrtrwt}}{numeric - square root of the weights applied
      to the residuals (read only)}
    \item{\code{variance}}{numeric - vector of unscaled variances (read
      only)}
    \item{\code{weights}}{numeric - prior weights (read only)}
    \item{\code{wrkResids}}{numeric - working residuals, on the scale of
      the linear predictor (read only)}
    \item{\code{wrkResp}}{numeric - working response, on the scale of
      the linear predictor (read only)}
    \item{\code{wrss}}{numeric scalar - weighted sum of squared 
      residuals (read only)}
    \item{\code{wtres}}{numeric - weighted residuals (read only)}
    \item{\code{y}}{numeric - response vector (read only)}
  }
}
\section{Methods}{
  \describe{
    \item{\code{updateMu(eta)}}{Install a new value of \code{eta} and
      update \code{mu} and derived quantities.  Note: the
      \code{sqrtrwt} vector is not updated but the \code{sqrtXwt} vector
      is.}
    \item{\code{updateWts()}}{Update \code{sqrtrwt} and \code{sqrtXwt}
      from the new value of the \code{variance}.}
  }
}

\details{
  After an object of this class has been constructed, it can be modified
  in two ways: the \code{updateMu} method takes a new value of the
  linear predictor, \code{eta}, and updates the mean, \code{mu}, and the
  weighted residuals,  \code{wtres}, returning the weighted residual sum
  of squares.  The \code{updateWts} method, which has no arguments,
  updates the variance and the residual and X weights (\code{sqrtrwt}
  and \code{sqrtXwt}.

  The implementation of penalized iteratively reweighted least squares
  (PIRLS) in \code{\link{glmer}} uses a Gauss-Newton algorithm to update
  the random-effects and optionally the fixed-effects parameters.  In
  this algorithm the weighted residuals, \code{wtres}, and the X
  weights, \code{sqrtXwts}, are all that is needed by the fixed-effects
  and random-effects modules to create a new linear predictor.  However,
  to determine an initial value for these coefficients we use one Fisher
  scoring iteration which requires the working response, \code{wrkResp},
  and the working weights, \code{sqrtWrkWt}.
}
%\source{}
%\references{}
\examples{
bb <- binomial()
set.seed(1)
eta0 <- sort(rnorm(10))
mu0 <- bb$linkinv(eta0)
y <- rbinom(length(mu0), 1, mu0)
gg <- new(glmerResp, bb, y)
gg$family
gg$link
gg$updateMu(eta0)  # returns weighted residual sum of squares
gg$updateWts()     # returns weighted residual sum of squares
comp <- c("y", "eta", "mu", "weights", "offset", "variance", "sqrtrwt",
          "wtres", "wrkResids", "wrkResp", "sqrtWrkWt")
names(comp) <- comp
sapply(comp, function(nm) gg[[nm]])
gg$residDeviance
}
\keyword{datasets}
