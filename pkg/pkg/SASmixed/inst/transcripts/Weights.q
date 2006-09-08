### $Id: Weights.q,v 1.1 1999/10/13 00:50:09 saikat Exp $
### Analysis of the weight-lifting program data given as data set 3.2(a)
### in "SAS System for Mixed Models"
options(contrasts = c(factor = "contr.SAS", ordered = "contr.poly"))
plot(Weights, layout = c(21,3), skip = rep(c(F,T,F,T,F),c(20,1,16,5,21)))
fm1Weight <- 
   lme( strength ~ Program * Time, data = Weights, random = ~ 1 | Subj,
       method = "ML" )
summary( fm1Weight )
summary( update( fm1Weight, method = "REML" ) )  # compare with output 3.1, p. 91
c( 3.0991, 1.0897 )^2
fm2Weight <- update( fm1Weight, random = ~ Time | Subj )
anova( fm1Weight, fm2Weight )
plot(augPred( fm2Weight ), layout = c(21,3),
      skip = rep(c(F,T,F,T,F),c(20,1,16,5,21)))
summary( fm2Weight )
fm3Weight <- update( fm2Weight, correlation = corAR1())
anova( fm2Weight, fm3Weight )
fm4Weight <- update( fm3Weight, strength ~ Program * (Time + I(Time^2)),
                    random = ~Time|Subj)
anova( fm1Weight, fm2Weight, fm3Weight, fm4Weight )
summary( fm4Weight )
intervals( fm4Weight )

