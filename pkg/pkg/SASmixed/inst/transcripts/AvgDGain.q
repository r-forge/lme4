### $Id: AvgDGain.q,v 1.1 1999/10/13 00:50:09 saikat Exp $
### Analysis of average daily weight gain data given as data set 5.3
### in "SAS System for Mixed Models"
options( contrasts = c(factor = "contr.SAS", ordered = "contr.poly") )
plot(AvgDailyGain)     # plot of adg versus Treatment by Block
fm1Adg <- lme( adg ~ InitWt * Treatment - 1, 
      data = AvgDailyGain, random = ~ 1 | Block, method = "ML" )
summary( fm1Adg ) # compare with output 5.1, p. 178
fm2Adg <- update( fm1Adg, adg ~ InitWt + Treatment )  # common slope model
anova( fm1Adg, fm2Adg )
summary( update( fm1Adg, adg ~ InitWt + Treatment - 1 ) )
