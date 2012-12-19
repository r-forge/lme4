load("randcrabdata.RData")
## randdata0: simulated data, in form suitable for plotting
## randdata: simulated data, in form suitable for analysis

fr  ## alive/dead formula
fr2 ## proportion alive formula (use with weights=initial.snail.density)

library(ggplot2)
library(grid)
zmargin <- theme(panel.margin=unit(0,"lines"))
theme_set(theme_bw())
g1 <- ggplot(randdata0,aes(x=snail.size,y=surv,colour=snail.size,fill=snail.size))+
    geom_hline(yintercept=1,colour="black")+
    stat_sum(aes(size=factor(..n..)),alpha=0.6)+
    facet_grid(.~ttt)+zmargin+
    geom_boxplot(fill=NA,outlier.colour=NULL,outlier.shape=3)+  ## set outliers to same colour as points
    ## (hard to see which are outliers, but it doesn't really matter in this case)
    scale_size_discrete("# obs",range=c(2,5))


library(lme4)
## FIXME: this test doesn't quite work, CRAN-lme4 is now 0.999999-0
if (packageVersion("lme4")>"0.999375-42") {
    ## using development lme4 ...
    try(glmer1 <- glmer(fr2,weights=initial.snail.density,family ="binomial", data=randdata))
    ## pwrssUpdate did not converge
    try(glmer1B <- glmer(fr,family ="binomial", data=randdata))
    if (require("lme4.0")) {
        detach("package:lme4",unload=TRUE)
        ## prop/weights formulation
        glmer1 <- glmer(fr2,weights=initial.snail.density,family ="binomial", data=randdata)
        ## alive/dead formulation
        glmer1B <- glmer(fr,family ="binomial", data=randdata)
        coef(glmer1B)
        fixef(glmer1B)
        detach("package:lme4.0")
    }
    if (require("glmmADMB")) {
        ## prop/weights formulation
        glmer1C <- glmmadmb(fr,family ="binomial", data=randdata)
    }
     if (require("MCMCglmm")) {
        ## prop/weights formulation
         ff <- lme4:::nobars(fr)
         ff <- ff[-2] ## delete response
         X <- model.matrix(ff,data=randdata)
         npar <- ncol(X)
        fit1D <- MCMCglmm(cbind(final.snail.density, snails.lost) ~
                          crab.speciesS + crab.speciesW + 
                          crab.sizeS + crab.sizeM + snail.sizeS + crab.speciesS:crab.sizeS + 
                          crab.speciesS:crab.sizeM + crab.speciesS:snail.sizeS +
                          crab.speciesW:snail.sizeS + 
                          crab.sizeS:snail.sizeS + crab.sizeM:snail.sizeS +
                          crab.speciesS:crab.sizeS:snail.sizeS +
                          crab.speciesS:crab.sizeM:snail.sizeS,
                  random=~us(1+snail.size):plot,
                          family="multinomial2",
                          data=randdata,
                          prior=list(B=list(mu=rep(0,npar),V=diag(npar)*1e3),
                                     G=list(list(nu=10,V=diag(2))),
                                     R=list(nu=10,V=1)),
                          verbose=FALSE)
    }
    if (require("coefplot2")) {
        cc2 <- coeftab(fit1D)
        rownames(cc2) <- gsub("^Sol\\.","",rownames(cc2))
        coefplot2(list(glmer1B,glmer1C,cc2),col=c(1,2,4))
    }
} else {
    ## CRAN version of lme4
    glmer1 <- glmer(fr2,weights=initial.snail.density,family ="binomial", data=randdata)
    glmer1B <- glmer(fr,family ="binomial", data=randdata)
}
    
fixef(glmer1)
fixef(glmer1B)
## note in this case that variances are not converging to zero
VarCorr(glmer1)

if (FALSE) {
    library(glmmADMB)
