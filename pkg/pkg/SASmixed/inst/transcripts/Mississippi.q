### $Id: Mississippi.q,v 1.1 1999/10/13 00:50:09 saikat Exp $
### Analysis of the Mississippi nitrogren concentrations given as data set
### 4.2 in "SAS System for Mixed Models"
options(contrasts=c(factor="contr.SAS", ordered="contr.poly"))
formula( Mississippi )
plot( Mississippi )
fm1Miss <- lme( y ~ 1, data = Mississippi, random = ~ 1 | influent,
               method = "ML")
summary( fm1Miss )         # compare with output 4.2, p. 143
fm1RMiss <- update( fm1Miss, method = "REML" )
summary( fm1RMiss )        # compare with output 4.1, p. 142
random.effects( fm1Miss )    # BLUP's of random effects on p. 144
random.effects( fm1Miss , aug = TRUE )   # including covariates
plot( random.effects( fm1Miss , aug = TRUE ), form = ~ Type )
random.effects( fm1RMiss )   # BLUP's of random effects on p. 142
intervals( fm1RMiss )        # interval estimates of variance components
c(2.9568, 7.9576, 21.416)^2 # compare to output 4.7, p. 148
fm2RMiss <- lme( y ~ Type, data = Mississippi, random = ~ 1 | influent,
      method = "REML" )
summary( fm2RMiss )         # compare to output 4.8 and 4.9, pp. 150-152
fm2Miss <- update( fm2RMiss, method = "ML" )  # get ML results too
anova( fm1Miss, fm2Miss )         # getting a p-value for the Type
  ## Notice that the p-value is considerably smaller than for the F-test.
