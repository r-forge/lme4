### $Id: Cultivar.q,v 1.1 1999/10/13 00:50:09 saikat Exp $
### Split-plot experiment with whole plots in randomized blocks (sec. 2.5)
options( contrasts = c(factor = "contr.SAS", ordered = "contr.poly") )
names( Cultivation )
formula( Cultivation )
plot( Cultivation, inner = ~ Inoc )
fm1Cult <- lme( drywt ~ Inoc * Cult, 
                 data = Cultivation, method = "ML",
                 random = list( Block = ~ 1, Cult = ~ 1 ) )
summary( fm1Cult )
logLik( update( fm1Cult, method = "REML" ) )  # check the same model is being fit
fm2Cult <- update( fm1Cult, drywt ~ Inoc + Cult )
fm3Cult <- update( fm1Cult, drywt ~ Inoc )
fm4Cult <- update( fm1Cult, drywt ~ 1 )
anova( fm1Cult, fm2Cult, fm3Cult, fm4Cult )
fm5Cult <- update( fm1Cult, drywt ~ Cult )
anova( fm1Cult, fm2Cult, fm5Cult, fm4Cult )
### AIC, BIC, and likelihood ratio tests all prefer fm3Cult
summary( update( fm3Cult, method = "REML" ) )  # REML estimates
