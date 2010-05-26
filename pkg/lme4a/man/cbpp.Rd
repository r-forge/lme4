\name{cbpp}
\alias{cbpp}
\docType{data}
\title{Contagious bovine pleuropneumonia}
\description{
  Contagious bovine pleuropneumonia (CBPP) is a major disease of cattle
  in Africa, caused by a mycoplasma.  This dataset describes the
  serological incidence of CBPP in zebu cattle during a follow-up survey
  implemented in 15 commercial herds located in the Boji district
  of Ethiopia.  The goal of the survey was to study the within-herd spread
  of CBPP in newly infected herds. Blood samples were quarterly
  collected from all animals of these herds to determine their CBPP
  status.  These data were used to compute the serological incidence of
  CBPP (new cases occurring during a given time period).  Some data are
  missing (lost to follow-up).
}
\usage{data(cbpp)
}
\format{
  A data frame with 56 observations on the following 4 variables.
  \describe{
    \item{\code{herd}}{A factor identifying the herd (1 to 15).}
    \item{\code{incidence}}{The number of new serological cases for a
      given herd and time period.}
    \item{\code{size}}{A numeric vector describing herd size at the
      beginning of a given time period.}
    \item{\code{period}}{A factor with levels \code{1} to \code{4}.}
  }
}
\details{
  Serological status was determined using a competitive enzyme-linked
  immuno-sorbent assay (cELISA).
}
\source{
  Lesnoff, M., Laval, G., Bonnet, P., Abdicho, S., Workalemahu, A.,
  Kifle, D., Peyraud, A., Lancelot, R., Thiaucourt, F. (2004)
  Within-herd spread of contagious bovine pleuropneumonia in Ethiopian highlands.
  \emph{Preventive Veterinary Medicine} \bold{64}, 27--40.
}
\examples{
## response as a matrix
(m1. <- glmer1(cbind(incidence, size - incidence) ~ period + (1 | herd),
               cbpp, binomial, verbose = TRUE))
(m1a <- glmer2(cbind(incidence, size - incidence) ~ period + (1 | herd),
               cbpp, binomial, verbose = 1))
## response as a vector of probabilities and usage of argument "weights"
m2. <- glmer1(incidence / size ~ period + (1 | herd), weights = size,
              cbpp, binomial, verbose = TRUE)
m2a <- glmer2(incidence / size ~ period + (1 | herd), weights = size,
              cbpp, binomial, verbose = 1)
## Confirm that these are equivalent:
stopifnot(all.equal(fixef(m1.), fixef(m2.), tol = 1e-7),
          all.equal(ranef(m1.), ranef(m2.), tol = 1e-7))
\dontshow{%% TODO: test logLik() working
  for(m in c(m1., m1a, m2., m2a)) {
    cat("-------\n\nCall: ", format(lme4a:::getCall(m)),"\n")
    print(logLik(m)); cat("AIC:", AIC(m), "\n") ; cat("BIC:", BIC(m),"\n")
  }
stopifnot(all.equal(logLik(m1.), logLik(m2.), tol = 1e-7),
          all.equal(AIC(m1.),    AIC(m2.),    tol = 1e-7),
          all.equal(BIC(m1.),    BIC(m2.),    tol = 1e-7))
if(FALSE) ## things are wrong "currently" (see below):
stopifnot(all.equal(logLik(m1.), logLik(m1a), tol = 0.02),
          all.equal(AIC(m1.),    AIC(m1a),    tol = 0.02),
          all.equal(BIC(m1.),    BIC(m1a),    tol = 0.02))
}

%% FIXME: Why are they so much different that we need such a high 'tol'  ????
%% DB: the glmer1() results are somewhat doubtful
## currently (2010-05-26) things are *really* wrong:
if(FALSE)
stopifnot(all.equal(fixef(m1.), fixef(m1a), tol = 0.02),
          all.equal(ranef(m2.), ranef(m2a), tol = 0.02))
if(FALSE)
stopifnot(all.equal(fixef(m1.), fixef(m1a), tol = 0.02),%% << FIXME tol !!
          all.equal(ranef(m2.), ranef(m2a), tol = 0.02))
}
\keyword{datasets}
