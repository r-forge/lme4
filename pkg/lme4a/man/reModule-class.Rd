\name{reModule-class}
\Rdversion{1.1}
\docType{class}
\alias{reModule-class}
\alias{reTrms-class}
\title{Class "reModule"}
\description{
  The \code{"reModule"} class incorporates information about the
  random effects in a mixed models.  The \code{"reTrms"} subclass
  includes information about random-effects terms from a model
  formula.
}
\section{Objects from the Class}{
  Objects from these classes are usually created as part of
  \code{\linkS4class{merMod}} object returned by functions
  \code{\link{lmer}}, \code{\link{glmer}} or \code{\link{nlmer}}.
}
\section{Slots}{
  \describe{
    \item{\code{L}:}{Sparse Cholesky factor of class \code{"CHMfactor"}} 
    \item{\code{Lambda}:}{Sparse matrix (class \code{"dgCMatrix"})
      representation of the relative covariance factor of the random
      effects}
    \item{\code{Lind}:}{integer vector of indices by which \code{theta}
      generates \code{Lambda}.  Must have 
      \code{length(Lind)==length(Lambda@x)} and \code{all(Lind \%in\%
	1:length(theta))}.  The update operation is of the form
      \code{Lambda@x[] = theta[Lind]}.}
    \item{\code{Zt}:}{Sparse model matrix, \code{"dgCMatrix"}, for the
      random effects.}
    \item{\code{lower}:}{Lower bounds on the elements of \code{theta}}
    \item{\code{theta}:}{Covariance parameter vector.}
    \item{\code{u}:}{Orthogonal random effects vector defined so that
      \code{b = Lambda \%*\% u}.}
    \item{\code{flist}:}{list of grouping factors from the
      random-effects terms (\code{"reTrms"} only).  Its \code{assign}
      attribute matches terms with grouping factors.}
    \item{\code{cnms}:}{list of column names from random-effects terms
      (\code{"reTrms"} only).}
  }
}
\section{Methods}{
  \describe{
    \item{getL}{\code{signature(x = "reModule")}: ... }
  }
}
%\references{
%}
\seealso{\code{\link{getL}}}
\examples{
showClass("reModule")
showClass("reTrms")
}
\keyword{classes}
