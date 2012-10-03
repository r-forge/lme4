do.plots <- FALSE
data(coalition2,package="Zelig")

form <- duration ~ invest + fract + polar + numst2 + crisis + (1 | country) 

tmpf <- function(m) {
  list(fixef=fixef(m),ranef=c(ranef(m)[[1]]),
       vc=c(unlist(VarCorr(m))),LL=logLik(m))
}

if (do.plots) {
pairs(coalition2[,1:7],pch=".",gap=0)
## ggplot2(coalition2,aes(x=fract,y=duration,
##                       colour=fract,shape=invest))+
##  geom_point()+facet_grid(cut_number())
}

## are there any other packages that will do Gamma?
require(lme4)
m_lme4 <- glmer(form, data   = coalition2, family = Gamma(link = log))
rlist <- list(lme4=tmpf(m_lme4))
packageVersion('lme4')
detach("package:lme4")

require(lme4a)
m_lme4a <- glmer(form, data   = coalition2, family = Gamma(link = log))
rlist <- c(rlist,list(lme4a=tmpf(m_lme4a)))
detach("package:lme4a")

require(lme4Eigen)
m_lme4e <- glmer(form, data   = coalition2, family = Gamma(link = log),
                 control=list(rhobeg=0.2, rhoend=2e-7),
                 tolPwrss=1e-8)
rlist <- c(rlist,list(lme4e=tmpf(m_lme4e)))
detach("package:lme4Eigen")

library(MASS)
m_glmmPQL <- glmmPQL(duration ~ invest + 
                        fract + 
                        polar + 
                        numst2 +
                        crisis,
                random = ~1 | country,
             data   = coalition2, 
             family = Gamma(link = log))
tt <- tmpf(m_glmmPQL)
## n.b. variance components don't match up exactly
tt$vc <- sum(diag(matrix(as.numeric(VarCorr(m_glmmPQL)),nrow=2)))
rlist <- c(rlist,list(glmmPQL=tt))

## detach("package:MASS")

rbind(sapply(rlist,"[[","fixef"),
      logLik=sapply(rlist,"[[","LL"),
      herdvar=sapply(rlist,"[[","vc"))

if (FALSE) {
  ## not working yet ...
  library(glmmADMB)

  f3 <- glmmadmb(duration ~ invest + 
                 fract + 
               polar + 
               numst2 +
               crisis,
               random = ~1 | country,
               data   = coalition2, 
               family = "Gamma", link="log",
               extra.args="-ams 160000000",
               save.dir="~/R/misc",
               verbose=TRUE)
}

## ./glmmadmb -maxfn 500 -noinit -shess -ams 120000000
