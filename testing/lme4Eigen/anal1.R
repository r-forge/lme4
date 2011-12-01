## utility functions for collecting bits of results

Sapply <- function(x,fun,...) {
  r <- lapply(x,fun,...)
  r0 <- r[[1]]
  sapply(r,function(x) if (is.null(x) || is.null(x[[1]]) || is.na(x)) rep(NA,length(r0)) else x)
}

getf <- function(x,n) {
  lapply(x,"[[",n)
}

pfun <- function(x) {
  paste(names(x),x,sep="=",collapse=",")
}
  
load("bobyqa1.RData")

## assemble data frame from results
opts <- getf(results,"options")
opts <- data.frame(rhobeg=Sapply(opts,"[[","rhobeg"),
                   rhoend=Sapply(opts,"[[","rhoend"),
                   tolPwrss=Sapply(opts,"[[","tolPwrss"))
deviance <- Sapply(results,"[[","deviance")
KKT <- t(Sapply(results,"[[","KKT"))
time <- unlist(getf(results,"time"))
pars <- do.call(rbind,getf(results,"parameters"))
dd <- data.frame(opts,deviance,KKT,pars,time)
dd$devOK <- with(dd,!is.na(deviance) &  (deviance-min(deviance))<1e-2)
dd$OKval <- with(dd,ifelse(devOK & isTRUE(KKT1) & isTRUE(KKT2),
                           "good",
                           ifelse(devOK,"OK","bad")))

library(ggplot2)

zmargin <- opts(panel.margin=unit(0,"lines"))
g0 <- ggplot(dd,aes(x=factor(log10(rhobeg/2)),
              y=factor(log10(rhoend/2)),colour=OKval))+
  geom_point(aes(size=time),alpha=0.8)+
  facet_grid(.~tolPwrss,labeller=label_both)+theme_bw()+zmargin+
  scale_size_continuous(limits=c(0.01,20),to=c(0.5,6))+
  labs(x=expression(paste("Start tol ",(log[10](rho/2)))),
       y=expression(paste("End tol ",(log[10](rho/2)))))

pdf("bobyqa1_plot.pdf",height=3,width=8)
print(g0)
dev.off()

with(dd,table(unlist(KKT1),unlist(KKT2),devOK))


