### $Id: Semi2.q,v 1.1 1999/10/13 00:50:09 saikat Exp $
### Analysis of the Oxide layer thickness data given as data set 4.4 in
### "SAS System for Mixed Models"
options(contrasts = c(factor = "contr.SAS", ordered = "contr.poly"))
formula(Semi2)
plot(Semi2)
fm1Semi2 <- lme( Thickness ~ 1, data = Semi2,
   random = ~ 1 | Lot/Wafer, method = "ML" )
summary( fm1Semi2 )
fm1RSemi2 <- update( fm1Semi2, method = "REML" )
summary( fm1RSemi2 )       # compare with output 4.13, p. 156
fm2Semi2 <- update( fm1Semi2, Thickness ~ Source )
anova( fm1Semi2, fm2Semi2 )
## Again, the p-value is smaller than that for the F test.
fm2RSemi2 <- update( fm2Semi2, method = "REML" )
summary( fm2RSemi2 )       # compare with output 4.15, p. 159
fm3Semi2 <- update( fm2Semi2, 
   random = list(Lot = pdDiag( ~ Source - 1 ), Wafer = ~ 1 ) )
summary( fm3Semi2 )
fm3RSemi2 <- update( fm3Semi2, method = "REML" )
summary( fm3RSemi2 )       # compare with output 4.17, p. 163
anova( fm1Semi2, fm2Semi2, fm3Semi2 )
