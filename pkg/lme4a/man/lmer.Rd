\name{lmer}
\alias{lmer}
\alias{lmer1}
\alias{lmer2}
\alias{glmer}
\alias{glmer1}
\alias{glmer2}
\alias{nlmer}
\alias{nlmer2}
\concept{LMM}% << so it's found; ditto for these:
\concept{GLMM}
\concept{NLMM}
\title{Fit Mixed-Effects Models}
\description{
  Fit a linear mixed model or a generalized linear mixed model or a
  nonlinear mixed model.
}
\usage{
lmer(formula, data, REML = TRUE, sparseX = FALSE,
      control = list(), start = NULL, verbose = 0,
      doFit = TRUE,%% TODO ??  add choice of optimizer as well (and compDev = TRUE/FALSE)  ??
      subset, weights, na.action, offset,
      contrasts = NULL, \dots)

glmer(formula, data, family = gaussian, sparseX = FALSE,
      control = list(), start = NULL,
      verbose = 0, doFit = TRUE, subset, weights,
      na.action, offset, contrasts = NULL, nAGQ = 1,
      mustart, etastart, \dots)

nlmer(formula, data, family = gaussian, sparseX = FALSE,
      start = NULL, verbose = 0, nAGQ = 1,
      doFit = TRUE, subset, weights, na.action, mustart, etastart,
      contrasts = NULL, control = list(), \dots)

##---------------------------------------------------------------------
## Somewhat experimental Variants :
##  -  lmer() is currently identical to  lmer2()
##  - glmer() is currently identical to glmer2() :

lmer1(formula, data, REML = TRUE, sparseX = FALSE,
      control = list(), start = NULL, verbose = 0,
      doFit = TRUE, compDev = TRUE,
      optimizer = c("nlminb", "bobyqa", "optimize"),%% MM: add "optim"; it did not get stuck initially like nlminb
      subset, weights, na.action, offset,
      contrasts = NULL, \dots)
lmer2(formula, data, REML = TRUE, sparseX = FALSE,
      control = list(), start = NULL, verbose = 0,
      doFit = TRUE,%% TODO ??  add choice of optimizer as well (and compDev = TRUE/FALSE)  ??
      subset, weights, na.action, offset,
      contrasts = NULL, \dots)

glmer1(formula, data, family = gaussian, sparseX = FALSE,
       compDev = TRUE, control = list(), start = NULL,
       verbose = 0, doFit = TRUE, subset, weights,
       na.action, offset, contrasts = NULL, nAGQ = 1,
       mustart, etastart, \dots)
glmer2(formula, data, family = gaussian, sparseX = FALSE,
       control = list(), start = NULL,
       verbose = 0, doFit = TRUE, subset, weights,
       na.action, offset, contrasts = NULL, nAGQ = 1,
       mustart, etastart, \dots)

nlmer2(formula, data, family = gaussian, start = NULL, verbose = 0,
       nAGQ = 1, doFit = TRUE, subset, weights, na.action, mustart,
       etastart, sparseX = FALSE, contrasts = NULL, control = list(),
       \dots)
}
\arguments{
  \item{formula}{a two-sided linear formula object describing the
    fixed-effects part of the model, with the response on the left of a
    \code{~} operator and the terms, separated by \code{+} operators, on
    the right.  The vertical bar character \code{"|"} separates an
    expression for a model matrix and a grouping factor.}
  \item{data}{an optional data frame containing the variables named in
    \code{formula}.  By default the variables are taken from the
    environment from which \code{lmer} is called.}
  \item{sparseX}{logical - should a sparse model matrix be used for the
    fixed-effects terms?  Defaults to \code{FALSE}.}
  \item{family}{a GLM family, see \code{\link[stats]{glm}} and
    \code{\link[stats]{family}}.}
  \item{REML}{logical - Should the estimates
    be chosen to optimize the REML criterion (as opposed to the
    log-likelihood)?  Defaults to \code{TRUE}.}
  \item{nAGQ}{a positive integer - the number of points per axis for
    evaluating the adaptive Gauss-Hermite approximation to the
    log-likelihood.  This defaults to 1, corresponding to the Laplacian
    approximation.  Values greater than 1 produce greater accuracy in
    the evaluation of the log-likelihood at the expense of speed.}
  \item{control}{a list of control parameters. See below for details.}
%% Fix this description
  \item{start}{a named list of starting values for the parameters in the
    model. If the list is of the same form as the \code{ST} slot, it is
    becomes the starting values of the \code{ST} slot.  If the list
    contains components named \code{fixef} and/or \code{ST}, these are
    used as the starting values for those slots.  (Setting starting
    values for \code{fixef} has no effect for a linear mixed model
    because the fixed-effects parameters do not appear in the profiled
    deviance.)  In \code{lmer} and \code{glmer} a numeric
    \code{start} argument of the appropriate length is used as the
    starting value of the parameter vector that determines the \code{ST}
    slot. In \code{nlmer} a numeric \code{start} argument is used as the
    starting values of the \code{fixef} slot.}
  \item{doFit}{logical scalar. When \code{doFit = FALSE} the model is
    not fit but instead a structure with the model matrices for the
    random-effects terms is returned, so they can be modified for
    special model forms. When \code{doFit = TRUE}, the default, the
    model is fit immediately.}
  \item{compDev}{logical - should compiled code be used for the deviance
    evaluation during the optimization of the parameters estimates?}
  \item{optimizer}{character - the name of the optimizer to use in
    determining the parameter estimates.  The default is
    \code{\link{nlminb}} from the \pkg{stats} package.  In some cases
    the \code{bobyqa} optimizer from the \pkg{minqa} package has
    proven faster and/or more reliable than \code{\link{nlminb}}.  For
    one-dimensional problems (a single, simple, scalar random-effects
    term), the \code{\link{optimize}} function can be used.}
  \item{subset, weights, na.action, offset, contrasts}{further model
    specification arguments as in \code{\link[stats]{lm}}; see there for
    details.}
  \item{mustart, etastart}{optional starting values on the scale of the
    conditional mean or the unbounded predictor as in
    \code{\link[stats]{glm}}; see there for details.}
  \item{verbose}{integer scalar.  If \code{> 0} verbose output is
    generated during the optimization of the parameter estimates.  If
    \code{> 1} verbose output is generated during the individual PIRLS
    steps.}
  \item{\dots}{other potential arguments.  A \code{method} argument was
    used in earlier versions of the package. Its functionality has been
    replaced by the \code{REML} and \code{nAGQ} arguments.}
}
\details{

  The \code{lmer} and \code{glmer} functions are nearly interchangeable.
  If \code{lmer} is called with a non-default \code{family} argument the
  call is replaced by a call to \code{glmer} with the current arguments.
  If \code{glmer} is called with the default \code{family}, namely the
  \code{\link[stats]{gaussian}} family with the identity link, then the call is
  replaced by a call to \code{lmer} with the current arguments.  (They
  are described as \dQuote{nearly} interchangeable because the
  \code{REML} argument only applies to calls to \code{lmer} and the
  \code{nAGQ} argument only applies to calls to \code{glmer}.)

  Additional standard arguments to model-fitting functions can be passed
  to \code{lmer} or \code{glmer} or \code{nlmer}:
  \describe{
    \item{subset}{an optional expression indicating the subset of the rows of
      \code{data} that should be used in the fit. This can be a logical
      vector, or a numeric vector indicating which observation numbers are
      to be included, or a  character  vector of the row names to be
      included.  All observations are included by default.}
    \item{na.action}{a function that indicates what should happen when the
      data contain \code{NA}s.  The default action (\code{na.fail}) prints
      an error message and terminate if there are any incomplete
      observations.}
    \item{control}{a named list of control parameters for the estimation
      algorithm, specifying only the ones to be changed from their
      default values.  Hence defaults to an empty list.\cr
      Possible control options and their default values are:
      \describe{
	\item{\code{msVerbose}:}{a logical value passed as the
	  \code{trace} argument to \code{nlminb} (see documentation on
	  that function).  Default is \code{getOption("verbose")}.}
	\item{\code{maxIter}:}{a positive integer passed as the
	  \code{maxIter} argument to \code{nlminb} (see documentation on
	  that function).  Default is \code{300}.}
	\item{\code{maxFN}:}{a positive integer specifying the
	  maximum number of evaluations of the deviance function allowed
	  during the optimization. Default is \code{900}}.
      }
    }
    \item{model, x}{logicals.  If \code{TRUE} the corresponding
      components of the fit (the model frame, the model matrices)
      are returned.}
  }

  The \code{lmer2} name exists only for backwards compatibility.
  Calling this function simply produces an equivalent call to \code{lmer}.
}
\value{
  An object of class \code{"\linkS4class{mer}"}, for which many methods
  are available.  See there for details.
}
\seealso{
  The \code{\linkS4class{mer}} class, \code{\link[stats]{lm}}
}
\examples{
## linear mixed models
(fm1. <- lmer1(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm2. <- lmer1(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
anova(fm1., fm2.)

(fm1 <- lmer2(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm2 <- lmer2(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
anova(fm1, fm2)

## dense vs sparse X  --------------------------
## use more sensical example !!
fm3.d <- lmer(Yield ~ 1|Batch, Dyestuff2)
fm3.s <- lmer(Yield ~ 1|Batch, Dyestuff2, sparseX=TRUE)#-> warning
## check "equality"
stopifnot(all.equal( coef(fm3.d), coef(fm3.s), tol = 1e-14),
          all.equal(sigma(fm3.d),sigma(fm3.s), tol = 1e-14),
          TRUE)


## generalized linear mixed model
(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial))
## nonlinear mixed models
\dontrun{
(nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~
              Asym + xmid + scal + (Asym|Tree),
              Orange,
              start = c(Asym = 200, xmid = 725, scal = 350)))
}
}
\keyword{models}
\keyword{methods}
