source("anal1.R")
## dd <- getdat("bobyqa1.RData")
## dd <- getdat("bobyqa_cbpp.RData")
## dd <- getdat("bobyqa_coalition2.RData")

dd <- getdat("Contraception_bobyqa1.RData")
dd2 <- getdat("Contraception_NM1.RData")

pdf("NM_tolpwrss_vs_deviance.pdf")

library(plyr)
library(reshape)
library(ggplot2)
dd2$logtolpwrss <- log10(dd2$opts2)
dd3 <- ddply(dd2,.(nAGQ),
      function(x) {
        transform(x,scdeviance=(deviance-min(deviance))/diff(range(deviance)))
      })
mm <- melt(subset(dd3,select=c(scdeviance,
                        deviance,logtolpwrss,time,nAGQ)),
           id.var=c("logtolpwrss","nAGQ"))
gg1 <- ggplot(mm,aes(x=logtolpwrss,y=value,colour=factor(nAGQ)))+
  geom_line(aes(group=nAGQ))+
              facet_grid(variable~.,scale="free")+theme_bw()+
  opts(panel.margin=unit(0,"lines"))
ggsave("lme4Eigen_NM_out1.pdf",gg1,device=pdf,width=5,height=6)

zmargin <- opts(panel.margin=unit(0,"lines"))
g0 <- ggplot(dd,aes(x=factor(log10(rhobeg/2)),
              y=factor(log10(rhoend/2)),colour=deviance))+
  geom_point(aes(size=time),alpha=0.8)+
  facet_grid(.~tolPwrss,labeller=label_both)+theme_bw()+zmargin+
  ## scale_size_continuous(limits=c(0.01,20),to=c(0.5,6))+
  labs(x=expression(paste("Start tol ",(log[10](rho/2)))),
       y=expression(paste("End tol ",(log[10](rho/2)))))

g0 %+% subset(dd,deviance<2384)
pdf("Contraception_bobyqa1_plot.pdf",height=3,width=8)
print(g0)
dev.off()

with(dd,table(unlist(KKT.1),unlist(KKT.2),devOK))

## c(gm1@ST[[1]],fixef(gm1)[1])
lme4sol <- data.frame(X1=0.64226046045739182411,
                      X2=-1.3985350510257590351,
                      deviance=100.09586074358657015)
## fixef(m_glmmADMB)[1], sqrt((c(m_glmmADMB$S[[1]])))
glmmADMBsol <- data.frame(X1=0.6443291,X2=-1.3985)

## where is the lme4/lme4a solution??
g2 <- ggplot(dd,aes(x=X1,y=X2,colour=deviance,shape=devOK))+
  geom_point(alpha=0.5)+theme_bw()

g2 %+% subset(dd,deviance<2385)
g2 +  geom_point(data=lme4sol,colour="green",shape=2)+
  geom_point(data=glmmADMBsol,colour="cyan",shape=2)
  

