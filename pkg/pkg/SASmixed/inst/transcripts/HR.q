### $Id: HR.q,v 1.1 1999/10/13 00:50:09 saikat Exp $
### Analysis of the Heart rate data given as data set 3.5 in
### "SAS System for Mixed Models"
options(contrasts = c(factor = "contr.SAS", ordered = "contr.poly"))
formula(HR)
plot(HR)                      # basic trellis plot
plot(HR, outer = ~ Drug )     # by drug type
fm1HR <- lme( HR ~ Time * Drug + baseHR, data = HR,  # linear trend in time
   random = ~ Time | Patient, method = "ML")
summary( fm1HR )
fm2HR <- update( fm1HR, weights = varPower(0) ) # use power-of-mean variance
summary( fm2HR )
intervals( fm2HR )             # variance function does not seem significant
anova( fm1HR, fm2HR )         # confirm with likelihood ratio
fm3HR <- update( fm1HR, HR ~ Time + Drug + baseHR ) # remove interaction
anova( fm1HR, fm3HR )
summary( fm3HR )
fm4HR <- update( fm3HR, HR ~ Time + baseHR )  # remove Drug term
anova( fm1HR, fm3HR, fm4HR )
summary( fm4HR )
plot( augPred( fm4HR, length = 2 ) )  # compare predictions and data
