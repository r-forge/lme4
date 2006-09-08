S> ### $Id: AvgDGain.St,v 1.1 1999/10/13 00:50:09 saikat Exp $
S> ### Analysis of average daily weight gain data given as data set 5.3
S> ### in "SAS System for Mixed Models"
S> options( contrasts = c(factor = "contr.SAS", ordered = "contr.poly") )
S> plot(AvgDailyGain)     # plot of adg versus Treatment by Block
S> fm1Adg <- lme( adg ~ InitWt * Treatment - 1, 
+       data = AvgDailyGain, random = ~ 1 | Block, method = "ML" )
S> summary( fm1Adg ) # compare with output 5.1, p. 178
Linear mixed-effects model fit by maximum likelihood
 Data: AvgDailyGain 
     AIC    BIC  logLik 
  29.971 44.628 -4.9855

Random effects:
 Formula:  ~ 1 | Block
        (Intercept) Residual 
StdDev:     0.47793  0.18723

Fixed effects: adg ~ InitWt * Treatment - 1 
                    Value Std.Error DF t-value p-value 
           InitWt  0.0045   0.00203 17  2.2285  0.0396
       Treatment0  0.4683   0.69774 17  0.6712  0.5111
      Treatment10  1.4435   0.62724 17  2.3013  0.0343
      Treatment20  0.5031   0.54174 17  0.9287  0.3661
      Treatment30  0.1751   0.75975 17  0.2304  0.8205
 InitWtTreatment0 -0.0023   0.00271 17 -0.8469  0.4088
InitWtTreatment10 -0.0035   0.00245 17 -1.4195  0.1738
InitWtTreatment20 -0.0012   0.00242 17 -0.4996  0.6238
 Correlation: 
                  InitWt Trtmn0 Trtm10 Trtm20 Trtm30 IntWT0 
       Treatment0  0.050                                   
      Treatment10 -0.032  0.054                            
      Treatment20  0.034  0.097  0.351                     
      Treatment30 -0.961  0.024  0.111  0.059              
 InitWtTreatment0 -0.780 -0.635  0.046 -0.024  0.750       
InitWtTreatment10 -0.807 -0.021 -0.530 -0.178  0.776  0.617
InitWtTreatment20 -0.856 -0.040 -0.107 -0.506  0.823  0.666
                  InWT10 
       Treatment0       
      Treatment10       
      Treatment20       
      Treatment30       
 InitWtTreatment0       
InitWtTreatment10       
InitWtTreatment20  0.775

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3    Max 
 -2.1282 -0.52113 -0.044351 0.53017 1.8544

Number of Observations: 32
Number of Groups: 8 
S> fm2Adg <- update( fm1Adg, adg ~ InitWt + Treatment )  # common slope model
S> anova( fm1Adg, fm2Adg )
       Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm1Adg     1 10 29.971 44.628 -4.9855                       
fm2Adg     2  7 27.617 37.877 -6.8086 1 vs 2  3.6463  0.3023
S> summary( update( fm1Adg, adg ~ InitWt + Treatment - 1 ) )
Linear mixed-effects model fit by maximum likelihood
 Data: AvgDailyGain 
     AIC    BIC  logLik 
  27.617 37.877 -6.8086

Random effects:
 Formula:  ~ 1 | Block
        (Intercept) Residual 
StdDev:     0.45883   0.2044

Fixed effects: adg ~ InitWt + Treatment - 1 
              Value Std.Error DF t-value p-value 
     InitWt 0.00277   0.00083 20  3.3418  0.0032
 Treatment0 0.25248   0.37756 20  0.6687  0.5113
Treatment10 0.73614   0.39019 20  1.8866  0.0738
Treatment20 0.71648   0.38265 20  1.8724  0.0758
Treatment30 0.80428   0.35582 20  2.2604  0.0351
 Correlation: 
            InitWt Trtmn0 Trtmn1 Trtm20 
 Treatment0 -0.859                     
Treatment10 -0.869  0.958              
Treatment20 -0.863  0.957  0.958       
Treatment30 -0.839  0.953  0.954  0.954

Standardized Within-Group Residuals:
     Min       Q1     Med      Q3    Max 
 -1.6326 -0.65331 0.09208 0.75471 2.0365

Number of Observations: 32
Number of Groups: 8 








