### $Id: Semiconductor.q,v 1.1 1999/10/13 00:50:09 saikat Exp $
### Analysis of the Semiconductor split-plot experiment, section 2.3 of
### "SAS System for Mixed Models"
options(contrasts = c(factor = "contr.SAS", ordered = "contr.poly"))
formula( Semiconductor )
plot( Semiconductor )
names( Semiconductor )
plot( Semiconductor, inner = ~ position )
fm1Semi <- lme( resistance ~ ET * position,
                 data = Semiconductor,
                 random = ~ 1 | Grp, method = "ML" )
summary( fm1Semi )
## check significance of the interaction
fm2Semi <- update( fm1Semi, resistance ~ ET + position )
fm3Semi <- update( fm1Semi, resistance ~ ET )
fm4Semi <- update( fm1Semi, resistance ~ position )
fm5Semi <- update( fm1Semi, resistance ~ 1 )
anova( fm5Semi, fm4Semi, fm2Semi, fm1Semi )
anova( fm5Semi, fm3Semi, fm2Semi, fm1Semi )
## AIC favors resistance ~ ET + position
## BIC favors resistance ~ 1
## Likelihood ratio seems to favor resistance ~ position

