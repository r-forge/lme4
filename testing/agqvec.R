data(cbpp,package="lme4")
getinfo.merMod <- function(m) c(fixef(m),deviance=deviance(m))

agqfun <- function(n=1,f=getinfo.merMod) {
    f(glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                         family = binomial, data = cbpp, nAGQ=n))
}

getpkg <- function(pkg) {
    paste(pkg,packageVersion(pkg),sep="_")
}
agqvec <- 1:25
dmat <- list()
pkg <- "lme4"
library(pkg,character.only=TRUE)
m <- t(sapply(agqvec,agqfun))
dmat[[1]] <- data.frame(pkg=getpkg(pkg),m,n=agqvec)
detach(paste0("package:",pkg),character.only=TRUE,unload=TRUE)

library(pkg,character.only=TRUE,lib.loc=".")
dmat[[2]] <- data.frame(pkg=getpkg(pkg),t(sapply(agqvec,agqfun)),n=agqvec)
detach(paste0("package:",pkg),character.only=TRUE,unload=TRUE)

pkg <- "lme4.0"
library(pkg,character.only=TRUE)
dmat[[3]] <- data.frame(pkg=getpkg(pkg),t(sapply(agqvec,agqfun)),n=agqvec)
detach(paste0("package:",pkg),character.only=TRUE,unload=TRUE)

pkg <- "glmmML"
library(pkg,character.only=TRUE)
getinfo.glmmML <- function(m) c(coef(m),deviance=deviance(m))
agqfun2 <- function(n=1) {
    getinfo.glmmML(glmmML(cbind(incidence, size - incidence) ~ period,
                   family = binomial,
                    cluster=herd,
                    method="ghq",
                    n.points=n,
                    data = cbpp))
}
dmat[[4]] <- data.frame(pkg=getpkg(pkg),t(sapply(agqvec,agqfun2)),n=agqvec)
detach(paste0("package:",pkg),character.only=TRUE,unload=TRUE)

dmat_all <- do.call(rbind,dmat)

library(ggplot2)
library(reshape2)
dmat2 <- melt(dmat_all,id.var=c("pkg","n"))
theme_update(theme_bw())
ggplot(dmat2,aes(x=n,y=value,colour=pkg))+facet_wrap(~variable,scale="free")+geom_line()+
    geom_point(aes(shape=pkg))
ggsave("agqtest1.pdf",width=10)

## if we want to plot deviations
aa1 <- acast(subset(dmat2,pkg!="lme4_0.999375.42"),pkg~n~variable)
aa2 <- sweep(aa1[1:2,,],c(2,3),aa1[3,,],"-")
dmat3 <- transform(setNames(melt(aa2),c("pkg","n","variable","value")),
                   variable=factor(variable,levels=levels(dmat2$variable)),
                   pkg=factor(pkg,levels=levels(dmat2$pkg)))
ggplot(dmat3,aes(x=n,y=value,colour=pkg))+facet_wrap(~variable,scale="free")+geom_line()+
    geom_point(aes(shape=pkg))+geom_hline(yintercept=0,colour="gray")
ggsave("agqtest2.pdf",width=10)
