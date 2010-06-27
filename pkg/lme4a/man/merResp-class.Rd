\name{merResp-class}
\Rdversion{1.1}
\docType{class}
\alias{merResp-class}
\alias{lmerResp-class}
\alias{glmerResp-class}
\alias{nlmerResp-class}
\alias{nglmerResp-class}
\title{"merResp" and derived classes}
\description{
  The \code{"merResp"} class is the virtual base class of
  response modules for mixed-effects models.  Classes that inherit from
  \code{"merResp"} include \code{\linkS4class{lmerResp}}, for linear
  mixed models, \code{\linkS4class{glmerResp}}, for generalized linear
  mixed models, and \code{\linkS4class{nlmerResp}}, for nonlinear mixed models.
}
\section{Objects from the Class}{
  Objects from these classes are usually created as part of
  \code{\linkS4class{merMod}} object returned by functions
  \code{\link{lmer}}, \code{\link{glmer}} or \code{\link{nlmer}}.
}
\section{Slots}{
  \describe{
    \item{\code{mu}:}{Fitted mean response.}
    \item{\code{offset}:}{offset in the linear predictor -- always
      present even if it is a vector of zeros.  In an
      \code{\linkS4class{nlmerResp}} object the length of the offset can be
      a multiple of the length of the response.}
    \item{\code{sqrtXwt}:}{the matrix of
    weights for the model matrices, derived from the \code{sqrtrwt} slot.}
    \item{\code{sqrtrwt}:}{Numeric vector of the square roots
      of the weights for the residuals.  For
      \code{\linkS4class{lmerResp}} and \code{\linkS4class{nlmerResp}}
      objects these are constant.  For \code{\linkS4class{glmerResp}}
      object these are updated at each iteration of the iteratively reweighted least
      squares algorithm.}
    \item{\code{weights}:}{Prior
      weights -- always present even when it is a vector of ones.}
    \item{\code{y}:}{Numeric response vector.}
    \item{\code{REML}:}{Logical indicator of whether to estimate
      parameters according to the REML criterion - \code{lmerResp}
      objects only.}
    \item{\code{family}:}{a glm family, see \code{\link{family}} for
      details - \code{glmerResp} objects only.}
    \item{\code{eta}:}{numeric vector, the linear predictor that is
      transformed to the conditional mean via the link function -
      \code{glmerResp} objects only.}
    \item{\code{n}:}{a numeric vector used for calculation of the aic family
      function (it is really only used with the binomial family but
      we need to include it everywhere) - \code{glmerResp} objects only.}
    \item{\code{nlenv}:}{an environment in which to evaluate the
      nonlinear model function - \code{nlmerResp} objects only.}
    \item{\code{nlmod}:}{an unevaluated call to the
      nonlinear model function - \code{nlmerResp} objects only.}
    \item{\code{pnames}:}{a character vector of parameter names -
      \code{nlmerResp} objects only.}
  }
}
\section{Methods}{
No methods defined with class "merResp" in the signature.
}
\examples{
showClass("merResp")
}
\keyword{classes}
