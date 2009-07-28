### $Id: CultMiss.q,v 1.1 1999/10/13 00:50:09 saikat Exp $
### A blocked split-plot with missing data (sec 2.7, pp. 68-75)
## Remove all observations in block 1, cultivar 'A'.
options( contrasts = c(factor = "contr.SAS", ordered = "contr.poly") )
CultMiss <- Cultivation[ Cultivation$Block != 1 | Cultivation$Cult != 'a', ]
dim(CultMiss)
plot( CultMiss, inner = ~ Inoc )
fm1CultM <- lme( drywt ~ Cult * Inoc, CultMiss, list(Block = ~ 1, Cult = ~ 1),
                method = "ML")
summary( fm1CultM )
fm2CultM <- update( fm1CultM, drywt ~ Cult + Inoc )
fm3CultM <- update( fm1CultM, drywt ~ Inoc )
fm4CultM <- update( fm1CultM, drywt ~ 1 )
anova( fm1Cult, fm2Cult, fm3Cult, fm4Cult )
# Essentially the same conclusions as for the balanced data
fm3RCultM <- update( fm3CultM, method = "REML" )
summary( fm3RCultM )
