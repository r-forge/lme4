### $Id: Multilocation.q,v 1.1 1999/10/13 00:50:09 saikat Exp $
### Analysis of the Multilocation data with fixed effects for the locations
options( contrasts = c(factor = "contr.SAS", ordered = "contr.poly") )
formula( Multilocation )
names( Multilocation )
### Create a Block %in% Location factor
Multilocation$Grp <-
  getGroups( Multilocation, form = ~ Location/Block, level = 2 )
fm1Mult <- lme( Adj ~ Location * Trt, data = Multilocation, ~ 1 | Grp,
               method = "ML")
summary( fm1Mult )
fm2Mult <- update( fm1Mult, Adj ~ Location + Trt )
fm3Mult <- update( fm1Mult, Adj ~ Location )
fm4Mult <- update( fm1Mult, Adj ~ Trt )
fm5Mult <- update( fm1Mult, Adj ~ 1 )
anova( fm1Mult, fm2Mult, fm3Mult, fm5Mult )
anova( fm1Mult, fm2Mult, fm4Mult, fm5Mult )
### AIC, BIC, and likelihood ratio tests all prefer model fm2Mult
summary( fm2Mult )
fm2RMult <- update( fm2Mult, method = "REML" ) # get REML estimates
summary( fm2RMult )
### Treating the location as a random effect
fm1MultR <- lme( Adj ~ Trt, data = Multilocation, method = "ML",
  random = list( Location = pdCompSymm( ~ Trt - 1 ), Block = ~ 1 ) )
summary( fm1MultR )
fm2MultR <- update( fm1MultR, random = list( Location = ~ Trt - 1, Block = ~ 1 ))
anova( fm1MultR, fm2MultR )
## No indication that a general variance-covariance is preferred to
## a compound symmetry structure.
fm1RMultR <- update( fm1MultR, method = "REML" )
summary( fm1RMultR )
c( 0.34116, 0.07497, 0.18596)^2  # compare with estimates, p. 84
