### $Id: PBIB.q,v 1.1 1999/10/13 00:50:09 saikat Exp $
### Analysis of the partially balanced incomplete blocked PBIB data
options( contrasts = c(factor = "contr.SAS", ordered = "contr.poly") )
formula( PBIB )
names( PBIB )
sapply( PBIB, data.class )
fm1PBIB <- lme( response ~ Treatment, data = PBIB, random = ~ 1 | Block,
               method = "ML")
summary( fm1PBIB )
plot( fm1PBIB, resid(.) ~ fitted(.) | Block )
plot( fm1PBIB, resid(.) ~ fitted(.) | Treatment, inner = ~ Block )
fm1RPBIB <- update( fm1PBIB, method = "REML" )
summary( fm1RPBIB )    # compare with output 1.7  pp. 24-25
## Testing for significant fixed effects for Treatment involved re-fitting
## the ML fit and comparing with anova().
fm2PBIB <- update( fm1PBIB, response ~ 1 )
anova( fm1PBIB, fm2PBIB )     
## The p-value is considerably smaller than that for the F test in PROC MIXED
