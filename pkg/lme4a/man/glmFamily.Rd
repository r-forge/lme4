\name{glmFamily}
\alias{glmFamily}
\docType{class}
\title{
 An C++ class for glm family objects
}
\description{
  An Rcpp module for a reference class representing a glm family
  object.  The constructor takes a \code{\link{family}} object.  The
  instance provides methods for functions in the family, making these
  available in \R{} and in \code{C++}.
}
\format{
  The format is:
Formal class 'C++Class' [package "Rcpp"] with 8 slots
  ..@ .Data       : chr "Rcpp_glmFamily"
  ..@ pointer     :<externalptr> 
  ..@ module      :<externalptr> 
  ..@ fields      :List of 2
  .. ..$ family:Formal class 'C++Field' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x6075748> 
  .. ..$ link  :Formal class 'C++Field' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x600b5c8> 
  ..@ methods     :List of 5
  .. ..$ devResid:Formal class 'C++OverloadedMethods' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x5fa9378> 
  .. ..$ linkFun :Formal class 'C++OverloadedMethods' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x5efe918> 
  .. ..$ linkInv :Formal class 'C++OverloadedMethods' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x5e67d38> 
  .. ..$ muEta   :Formal class 'C++OverloadedMethods' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x5dd4798> 
  .. ..$ variance:Formal class 'C++OverloadedMethods' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x5d5dfb0> 
  ..@ constructors:List of 1
  .. ..$ :Formal class 'C++Constructor' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x581e370> 
  ..@ generator   :Formal class 'refObjectGenerator' [package "methods"] with 1 slots
  .. .. ..@ .xData:<environment: 0x5d57620> 
  ..@ docstring   : chr ""
}
\details{
  A glm \code{\link{family}} object is a list of functions including the
  link function (\code{linkfun}), the inverse link function
  (\code{linkinv}), the \code{variance} function, the \code{mu.eta}
  function that calculates the derivative of the mean with respect to
  the linear predictor, and a function that evaluates the (squared)
  deviance residuals.

  A \code{Rcpp_glmFamily} object allows these functions to be accessed
  within \code{C++} code or within \R{} code.

  For several common families the transformation functions:
  \code{linkFun}, \code{linkInv}, \code{muEta} and \code{variance}, are
  evaluated in compiled code.  }

\note{
  This class is primarily used for testing.  It is not used directly.
}
\seealso{\code{\link{glmerResp}}, \code{\link{glmer}}}
%\source{}
%\references{}
\examples{
bFam <- new(glmFamily, binomial())
bFam$family
bFam$link
set.seed(1)
mu <- sort(runif(10))	        # mean vector with elements in (0,1)
eta <- bFam$linkFun(mu)         # linear predictor
cbind(mu, eta, muEta=bFam$muEta(eta), variance=bFam$variance(mu))
}
\keyword{methods}
