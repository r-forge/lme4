\name{profile-methods}
\docType{methods}
\alias{profile-methods}
\alias{profile,ANY-method}
\alias{profile,lmer-method}
\alias{profile,lmerenv-method}
\title{Methods for profile() of [ng]lmer Objects}
\description{
  Methods for function \code{\link{profile}} (package \pkg{stats}),
  here for profiling (fitted) mixed effect models.
}
\usage{
\S4method{profile}{lmerenv}(fitted, alphamax = 0.01, maxpts = 100,
          delta = cutoff/8, tr = 0, \dots)
}
\arguments{
  \item{fitted}{a fitted model, e.g., the result of \code{\link{lmer}(..)}.}
  \item{alphamax}{used when \code{delta} is unspecified, as probability
    ... to compute \code{delta} ...}
  \item{maxpts}{ ... % FIXME
}
  \item{delta}{ ... % FIXME
}
  \item{tr}{ ... % FIXME
}
  \item{\dots}{potential further arguments for \code{profile} methods.}
}

\section{Methods}{
  FIXME: These (signatures) will change soon --- document \bold{after} change!

  \describe{
    \item{\code{signature(fitted = "lmer")}}{ ...
    }

    \item{\code{signature(fitted = "lmerenv")}}{ ...
    }
  }
}
\examples{
%% Do keep at least *one* such example! -- this is (also) a regression test
fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)
tpr <- profile(fm01ML)
(confint(tpr) -> CIpr)
print(xyplot(tpr))
stopifnot(dim(CIpr) == c(3,2),
          all.equal(unname(CIpr[2,]),c(3.64362, 4.21446), tol=1e-6))
}
\keyword{methods}
