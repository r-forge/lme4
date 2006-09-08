### $Id: Bond.q,v 1.1 1999/10/13 00:50:09 saikat Exp $
## Set options 
options( prompt = "S> ", digits = 5, width = 65, 
  contrasts = c(factor = "contr.SAS", ordered = "contr.poly") )
plot(Bond)                       # dotplot by Ingot
plot(Bond, inner = ~ Metal)      # different symbols for different Metals
formula(Bond)                    # check the formula
fm1Bond <- lme( pressure ~ Metal, data = Bond, random = ~ 1 | Ingot,
               method = "ML")
summary( fm1Bond )
## default criterion in lme is maximum likelihood (ML).
## Re-fit to get REML results
fm1RBond <- update( fm1Bond, method = "REML" )
summary( fm1RBond )              # compare with output 1.1 on p. 6
logLik( fm1RBond )                # log-restricted-likelihood
c(3.3835,3.2205)^2                # variance estimates 
##  To test the need for the Metal term in the fixed effects, 
##  re-fit and use anova.  You must use ML to do this.  RML results are
##  not comparable if the fixed-effects specification changes.
fm2Bond <- update( fm1Bond, pressure ~ 1 )
anova( fm1Bond, fm2Bond )
