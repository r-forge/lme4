library(optimx)
library(lme4Eigen)
form <- cbind(incidence, size - incidence) ~ period + (1 | herd)

ff <- glmer(form, cbpp, binomial, devFunOnly=2L)
environment(ff)$tolPwrss                # current default (but should be changed)
optimx(with(environment(ff), c(pp$theta, pp$beta0)), ff,
       lower=c(0, rep.int(-Inf, 4L)), control=list(all.methods=TRUE))

ff <- glmer(form, cbpp, binomial, devFunOnly=2L, tolPwrss=1e-9)
optimx(with(environment(ff), c(pp$theta, pp$beta0)), ff,
       lower=c(0, rep.int(-Inf, 4L)), control=list(all.methods=TRUE))

ff <- glmer(form, cbpp, binomial, devFunOnly=2L, tolPwrss=1e-10)
optimx(with(environment(ff), c(pp$theta, pp$beta0)), ff,
       lower=c(0, rep.int(-Inf, 4L)), control=list(all.methods=TRUE))

