
## get data
if (!file.exists("search.csv")) {
    download.file("http://psych.colorado.edu/~westfaja/search.csv",
                  dest="search.csv")
}
                
search <- read.csv("search.csv")

## set contrasts
contrasts(search$target) <- contr.helmert(2)
contrasts(search$val) <- cbind(negVsPos=-1:1,
                               neuVsAll=c(1,-2,1))

form <- Intensity ~ target*val + (target*val|sub) + (target|image)
## fit model
library(lme4)
t1 <- system.time( mod1 <- lmer(form, data=search))

## try with more iterations
t2 <-system.time( mod2 <- lmer(Intensity ~ form,
                               data=search,
                               control=list(maxfun=100000)))

require(optimx)
t3 <- system.time( mod3 <- lmer(form, data=search,
                                optimizer="optimx",
                                control=list(method="nlminb")))             

t5 <- system.time( mod5 <- lmer(form, data=search,
                                optimizer="bobyqa"))


detach("package:lme4")
library(lme4.0)
t4 <- system.time( mod4 <- lmer(form, data=search))
detach("package:lme4.0")

library(bbmle)
AICtab(mod1,mod2,mod3,mod4,mod5)
sapply(list(t1,t2,t3,t4,t5),"[","elapsed")

library(lme4)
source("slice.R")
mm <- update(mod3,devFunOnly=TRUE)
tt <- getME(mod3,"theta")
s1 <- slice(tt,mm,dim=2,verbose=TRUE)

png("westfall_splom.png",height=1500,width=1500)
splom(s1)
dev.off()


library(coefplot2)  ## from r-forge
coeftab(mod1)
coeftab(mod2)
coeftab(mod3)
coeftab(mod4)

getvcov <- function(x) {
    cc <- suppressWarnings(coeftab(x,ptype="vcov"))
    cc1 <- cc[,1]
    names(cc1) <- rownames(cc)
    cc1
}

sapply(list(mod1,mod2,mod3,mod4,mod5),getvcov)
sapply(list(mod1,mod2,mod3,mod5),fixef)

savestuff <- c(ls(pattern="(t[0-9]|mod[0-9])"),"s1","mm","tt","search")
save(list=savestuff,
     file="westfall.RData")
