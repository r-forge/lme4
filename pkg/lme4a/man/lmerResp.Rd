\name{lmerResp}
\alias{lmerResp}
\docType{class}
\title{
  An C++ class for \code{lmerResp} objects
}
\description{
  The \code{\linkS4class{lmerResp}} S4 class represents the information
  associated with the response in an \code{\link{lm}} or
  \code{\link{lmer}} model. This C++Class encapsulates the information
  and provides methods accessible in \R{} or in \code{C++}.
}
\format{
  The format is:
Formal class 'C++Class' [package "Rcpp"] with 8 slots
  ..@ .Data       : chr "Rcpp_lmerResp"
  ..@ pointer     :<externalptr> 
  ..@ module      :<externalptr> 
  ..@ fields      :List of 8
  .. ..$ mu     :Formal class 'C++Field' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x34a5e18> 
  .. ..$ offset :Formal class 'C++Field' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x16509c8> 
  .. ..$ sqrtXwt:Formal class 'C++Field' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x4e84600> 
  .. ..$ sqrtrwt:Formal class 'C++Field' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x4ee22a8> 
  .. ..$ weights:Formal class 'C++Field' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x4f42780> 
  .. ..$ wrss   :Formal class 'C++Field' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x4fa3e58> 
  .. ..$ wtres  :Formal class 'C++Field' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x5004a60> 
  .. ..$ y      :Formal class 'C++Field' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x5063ab0> 
  ..@ methods     :List of 2
  .. ..$ updateMu :Formal class 'C++OverloadedMethods' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x50c0c88> 
  .. ..$ updateWts:Formal class 'C++OverloadedMethods' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x51448f8> 
  ..@ constructors:List of 4
  .. ..$ :Formal class 'C++Constructor' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x51cacd0> 
  .. ..$ :Formal class 'C++Constructor' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x522a2e8> 
  .. ..$ :Formal class 'C++Constructor' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x528b8e0> 
  .. ..$ :Formal class 'C++Constructor' [package "Rcpp"] with 1 slots
  .. .. .. ..@ .xData:<environment: 0x52ede58> 
  ..@ generator   :Formal class 'refObjectGenerator' [package "methods"] with 1 slots
  .. .. ..@ .xData:<environment: 0x52b5408> 
  ..@ docstring   : chr ""
}
\details{
  The methods in this class parallel those in the
  \code{"\link{glmerResp}"} class but typically \code{updateWts} is
  not called.  The \code{updateMu} method is called after every change
  in the variance component parameter, providing an updated value of the
  residual sum of squares.

  The first argument to the constructors (other than from an S4 object)
  is the REML indicator, which is an integer value of 0, for
  ML estimates, or \code{p}, the rank of the fixed-effects model matrix,
  for REML estimates.
}
%\source{}
%\references{}
\examples{
ll <- with(sleepstudy, new(lmerResp, 0L, Reaction))
comp <- c("y", "mu", "weights", "offset", "sqrtrwt", "wtres")
names(comp) <- comp
head(sapply(comp, function(nm) ll[[nm]]), n = 10)
ll$wrss
}
\keyword{methods}
