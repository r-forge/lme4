\name{feModule-class}
\Rdversion{1.1}
\docType{class}
\alias{feModule-class}
\alias{deFeMod-class}
\alias{spFeMod-class}
\title{Class "deFeMod"}
\description{
  The \code{"feModule"} class is the virtual base class of
  fixed-effects modules for mixed-effects models.  Classes that inherit from
  \code{"feModule"} are \code{"deFeMod"} and \code{"spFeMod"}, the
  distinction being whether the model matrix is a dense ("de") or a
  sparse ("sp") matrix.
}
\section{Objects from the Class}{
  Objects from these classes are usually created as part of
  \code{\linkS4class{merMod}} object returned by functions
  \code{\link{lmer}}, \code{\link{glmer}} or \code{\link{nlmer}}.
}
\section{Slots}{
  \describe{
    \item{\code{RZX}:}{The off-diagonal part of the Cholesky
      decomposition for the penalized least squares problem - of class 
      \code{"\linkS4class{dgeMatrix}"} for \code{"deFeMod"} objects or
      \code{"\linkS4class{dgCMatrix}"} for \code{"spFeMod"} objects.}
    \item{\code{RX}:}{The lower-right part of the Cholesky
      decomposition for the penalized least squares problem - of class 
      \code{"\linkS4class{Cholesky}"} for \code{"deFeMod"} objects or
      \code{"\linkS4class{CHMfactor}"} for \code{"spFeMod"} objects.}
    \item{\code{X}:}{Model matrix for the fixed-effects parameters - of
      class \code{"\linkS4class{dgeMatrix}"} for \code{"deFeMod"}
      objects or \code{"\linkS4class{dgCMatrix}"} for \code{"spFeMod"}
      objects.} 
    \item{\code{beta}:}{Numeric vector of fixed-effects parameters}
  }
}
\section{Methods}{
No methods defined with class "feModule" in the signature.
}
\seealso{\code{\linkS4class{merMod}}}
\examples{
showClass("deFeMod")
}
\keyword{classes}
