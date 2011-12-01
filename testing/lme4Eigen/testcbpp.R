library(lme4)
cbpp <- cbpp
form <- cbind(incidence, size - incidence) ~ period + (1 | herd)

## utility for extracting *mer info

tmpf <- function(m) {
  list(fixef=fixef(m),ranef=c(ranef(m)[[1]]),vc=c(unlist(VarCorr(m))),LL=logLik(m))
}

## 1. lme4
m_lme4 <- glmer(form, family = binomial, data = cbpp)
rlist <- list(lme4=tmpf(m_lme4))
detach("package:lme4")

## 2. lme4a
library(lme4a)
m_lme4a <- glmer(form, family = binomial, data = cbpp)
rlist <- c(rlist,list(lme4a=tmpf(m_lme4a)))
detach("package:lme4a")

## 3. lme4eigen
library(lme4Eigen)

## 3A. default settings
m_lme4e <- glmer(form, family = binomial, data = cbpp)
rlist <- c(rlist,list(lme4e=tmpf(m_lme4e)))

## 3B. change rhobeg/rhoend settings back to lme4a defaults
m_lme4eB <- glmer(form,family = binomial, data = cbpp,
                  control=list(rhobeg=0.2, rhoend=2e-7))
rlist <- c(rlist,list(lme4eB=tmpf(m_lme4eB)))

## 3C. rhobeg/rhoend settings plus decreased value of tolPwrss
m_lme4eC <- glmer(form,
            family = binomial, data = cbpp, control=list(rhobeg=0.2, rhoend=2e-7),
                  tolPwrss=1e-8)
rlist <- c(rlist,list(lme4eC=tmpf(m_lme4eC)))

detach("package:lme4Eigen")

## 4. glmmADMB
## n.b. glmmADMB needs to be installed from r-forge -- not on CRAN
library(glmmADMB)

m_glmmADMB <- glmmadmb(cbind(incidence, size - incidence) ~ period + (1 | herd),
            family = "binomial", data = cbpp)
rlist <- c(rlist,list(glmmADMB=tmpf(m_glmmADMB)))
detach("package:glmmADMB")

## 5. glmmML
library(glmmML)
m_glmmML <- glmmML(cbind(incidence, size - incidence) ~ period,
                   cluster=herd, family="binomial", data=cbpp)
rlist <- c(rlist,list(glmmML=list(fixef=coef(m_glmmML),ranef=m_glmmML$posterior.modes,
                        LL=-m_glmmML$deviance/2,vc=m_glmmML$sigma^2)))

rbind(sapply(rlist,"[[","fixef"),
      logLik=sapply(rlist,"[[","LL"),
      herdvar=sapply(rlist,"[[","vc"))
## sapply(rlist,"[[","ranef")



