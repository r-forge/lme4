\name{getBounds}
\Rdversion{1.1}
\alias{getBounds}
\alias{getBounds-methods}
\alias{getBounds,environment-method}
\alias{getBounds,optenv-method}
\alias{getPars}
\alias{getPars-methods}
\alias{getPars,environment-method}
\alias{getPars,optenv-method}
\alias{setPars}
\alias{setPars-methods}
\alias{setPars,environment,numeric-method}
\alias{setPars,optenv,numeric-method}
\alias{env}
\alias{env,optenv-method}
\title{
  Generic functions required by optimization environment objects
}
\description{
  Methods for these three generic functions are used to define an
  optimization environment class.  The state of the system is
  encapsulated in the environment and the optimization is performed by
  repeated calls to the \code{setPars} method which sets a new value of
  the parameter values and returns the value of the objective function
  at these parameter values.
}
\usage{
getBounds(x, \dots)
getPars(x, \dots)
setPars(x, pars, \dots)
env(x, \dots)
}
\arguments{
  \item{x}{an object of a class, such as \code{\linkS4class{optenv}},
    for which suitable methods are defined.}
  \item{pars}{numeric vector of parameter values to be installed in the
    object, to be followed by an evaluation of the objective function.
    The length of the parameter vector must match the length of the
    value returned by the \code{getPars} method and the numeric values
    should be a feasible point for the bounds returned by the
    \code{getBounds} method.}
  \item{\dots}{optional additional arguments (none used at present).}
}
\value{
  A \code{getPars} method returns the numeric vector of the object's
  parameters that are to be optimized.  A \code{getBounds} method returns
  a matrix with two columns, the lower and upper bounds on the parameter
  values.  A \code{setPars} method sets the value of the parameters in
  \code{x} to \code{pars} and returns the value of the objective function
  at these parameter values.

  An \code{env} method returns the common environment accessed by the other
  methods.
}
\section{Methods}{
\describe{
  \item{x = "environment"}{\code{getBounds} or \code{getPars} method}
  \item{x = "optenv"}{\code{getBounds} or \code{getPars} or \code{env} method}
  \item{x = "environment", pars = "numeric"}{\code{setPars} method}
  \item{x = "optenv", pars = "numeric"}{\code{setPars} method}
}}
\keyword{methods}
