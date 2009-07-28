### $Id: Demand.q,v 1.1 1999/10/13 00:50:09 saikat Exp $
### Analysis of the per capita demand deposits data given as data
### set 3.6 in "SAS System for Mixed Models"
names( Demand )
formula( Demand )
# the "grp" factor is a dummy factor with only one level.
unique( Demand$grp )
# Crossed random-effects factors have to be created by pdIdent applied to
# the indicator variables and joined together by pdBlocked.
fm1Demand <- lme( log(d) ~ log(y) + log(rd) + log(rt) + log(rs), data = Demand,
  random = list(grp = pdBlocked(list(pdIdent(~ State - 1), pdIdent(~ Year - 1)))),
  method = "REML" )
summary( fm1Demand )        # compare to output 3.13, p. 132
