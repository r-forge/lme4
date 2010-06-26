\name{env}
\Rdversion{1.1}
\alias{env}
\alias{env,optenv-method}
\title{
  Environment of an optimization environment
}
\description{
  The slots of an \code{\linkS4class{optenv}} object are functions that
  share a parent environment.  The \code{env} extractor accesses this
  environment.
}
\usage{
env(x, \dots)
}
\arguments{
  \item{x}{an object of a class, such as \code{\linkS4class{optenv}},
    for which suitable methods are defined.}
  \item{\dots}{optional additional arguments (none used at present).}
}
\value{
  An \code{env} method returns the common environment accessed by the other
  methods.
}
\section{Methods}{
\describe{
  \item{x = "optenv"}{\code{env} method}
}}
\keyword{methods}
