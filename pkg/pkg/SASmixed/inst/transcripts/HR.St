S> ### $Id: HR.St,v 1.1 1999/10/13 00:50:09 saikat Exp $
S> ### Analysis of the Heart rate data given as data set 3.5 in
S> ### "SAS System for Mixed Models"
S> options(contrasts = c(factor = "contr.SAS", ordered = "contr.poly"))
S> formula(HR)
HR ~ Time | Patient
S> plot(HR)                      # basic trellis plot
S> plot(HR, outer = ~ Drug )     # by drug type
S> fm1HR <- lme( HR ~ Time * Drug + baseHR, data = HR,  # linear trend in time
+    random = ~ Time | Patient, method = "ML")
S> summary( fm1HR )
Linear mixed-effects model fit by maximum likelihood
 Data: HR 
     AIC    BIC  logLik 
  809.75 840.42 -393.88

Random effects:
 Formula:  ~ Time | Patient
 Structure: General positive-definite
            StdDev   Corr 
(Intercept) 7.0549 (Inter
       Time 5.3200 -0.581
   Residual 4.9353       

Fixed effects: HR ~ Time * Drug + baseHR 
              Value Std.Error DF t-value p-value 
(Intercept)  33.976    9.6799 93  3.5099  0.0007
       Time  -3.197    2.9744 93 -1.0749  0.2852
      Druga   3.599    4.0114 20  0.8973  0.3803
      Drugb   7.091    3.9909 20  1.7768  0.0908
     baseHR   0.543    0.1093 20  4.9737  0.0001
  TimeDruga  -7.501    4.2064 93 -1.7833  0.0778
  TimeDrugb  -3.989    4.2064 93 -0.9484  0.3454
 Correlation: 
          (Intr)   Time  Druga  Drugb baseHR TimDrg 
     Time -0.166                                   
    Druga -0.309  0.400                            
    Drugb -0.245  0.402  0.501                     
   baseHR -0.957  0.000  0.109  0.041              
TimeDruga  0.117 -0.707 -0.566 -0.285  0.000       
TimeDrugb  0.117 -0.707 -0.283 -0.569  0.000  0.500

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3    Max 
 -2.1005 -0.48442 0.0059578 0.57901 2.3471

Number of Observations: 120
Number of Groups: 24 
S> fm2HR <- update( fm1HR, weights = varPower(0) ) # use power-of-mean variance
S> summary( fm2HR )
Linear mixed-effects model fit by maximum likelihood
 Data: HR 
     AIC    BIC  logLik 
  809.02 842.47 -392.51

Random effects:
 Formula:  ~ Time | Patient
 Structure: General positive-definite
              StdDev   Corr 
(Intercept) 7.097978 (Inter
       Time 5.769233 -0.553
   Residual 0.031682       

Variance function:
 Structure: Power of variance covariate
 Formula:  ~ fitted(.) 
 Parameter estimates:
  power 
 1.1471
Fixed effects: HR ~ Time * Drug + baseHR 
              Value Std.Error DF t-value p-value 
(Intercept)  34.122    9.6552 93  3.5341  0.0006
       Time  -3.818    2.9847 93 -1.2793  0.2040
      Druga   2.864    4.0128 20  0.7137  0.4837
      Drugb   6.646    4.0183 20  1.6540  0.1137
     baseHR   0.543    0.1091 20  4.9793  0.0001
  TimeDruga  -6.014    4.1681 93 -1.4429  0.1524
  TimeDrugb  -2.882    4.2829 93 -0.6728  0.5027
 Correlation: 
          (Intr)   Time  Druga  Drugb baseHR TimDrg 
     Time -0.174                                   
    Druga -0.312  0.397                            
    Drugb -0.246  0.395  0.498                     
   baseHR -0.956  0.011  0.111  0.043              
TimeDruga  0.121 -0.716 -0.563 -0.283 -0.004       
TimeDrugb  0.113 -0.697 -0.275 -0.565  0.001  0.499

Standardized Within-Group Residuals:
     Min       Q1       Med     Q3    Max 
 -2.1611 -0.48018 0.0035072 0.5435 2.1633

Number of Observations: 120
Number of Groups: 24 
S> intervals( fm2HR )             # variance function does not seem significant
Approximate 95% confidence intervals

 Fixed effects:
                lower      est.    upper 
(Intercept)  15.51646  34.12222 52.72799
       Time  -9.56996  -3.81839  1.93318
      Druga  -5.25875   2.86395 10.98666
      Drugb  -1.48758   6.64625 14.78008
     baseHR   0.32236   0.54317  0.76398
  TimeDruga -14.04611  -6.01403  2.01804
  TimeDrugb -11.13503  -2.88173  5.37157

 Random Effects:
  Level: Patient 
                         lower     est.     upper 
      sd((Intercept))  4.97762  7.09798 10.121554
             sd(Time)  3.20638  5.76923 10.380560
cor((Intercept),Time) -0.85589 -0.55318  0.031879

 Variance function:
          lower   est.  upper 
power -0.085702 1.1471 2.3799

 Within-group standard error:
      lower     est.  upper 
 0.00014209 0.031682 7.0642
S> anova( fm1HR, fm2HR )         # confirm with likelihood ratio
      Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm1HR     1 11 809.75 840.42 -393.88                       
fm2HR     2 12 809.02 842.47 -392.51 1 vs 2  2.7355  0.0981
S> fm3HR <- update( fm1HR, HR ~ Time + Drug + baseHR ) # remove interaction
S> anova( fm1HR, fm3HR )
      Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm1HR     1 11 809.75 840.42 -393.88                       
fm3HR     2  9 808.92 834.01 -395.46 1 vs 2  3.1651  0.2055
S> summary( fm3HR )
Linear mixed-effects model fit by maximum likelihood
 Data: HR 
     AIC    BIC  logLik 
  808.92 834.01 -395.46

Random effects:
 Formula:  ~ Time | Patient
 Structure: General positive-definite
            StdDev   Corr 
(Intercept) 7.2461 (Inter
       Time 6.1343 -0.604
   Residual 4.9359       

Fixed effects: HR ~ Time + Drug + baseHR 
              Value Std.Error DF t-value p-value 
(Intercept)  36.049    9.5130 95  3.7895  0.0003
       Time  -7.027    1.8176 95 -3.8663  0.0002
      Druga  -0.453    3.2773 20 -0.1383  0.8914
      Drugb   4.936    3.2526 20  1.5175  0.1448
     baseHR   0.543    0.1083 20  5.0171  0.0001
 Correlation: 
       (Intr)   Time  Druga  Drugb 
  Time -0.103                     
 Druga -0.297  0.000              
 Drugb -0.219  0.000  0.502       
baseHR -0.965  0.000  0.132  0.050

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3    Max 
 -2.1262 -0.54393 -0.031224 0.45798 2.4581

Number of Observations: 120
Number of Groups: 24 
S> fm4HR <- update( fm3HR, HR ~ Time + baseHR )  # remove Drug term
S> anova( fm1HR, fm3HR, fm4HR )
      Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm1HR     1 11 809.75 840.42 -393.88                       
fm3HR     2  9 808.92 834.01 -395.46 1 vs 2  3.1651  0.2055
fm4HR     3  7 808.21 827.72 -397.11 2 vs 3  3.2925  0.1928
S> summary( fm4HR )
Linear mixed-effects model fit by maximum likelihood
 Data: HR 
     AIC    BIC  logLik 
  808.21 827.72 -397.11

Random effects:
 Formula:  ~ Time | Patient
 Structure: General positive-definite
            StdDev   Corr 
(Intercept) 7.6072 (Inter
       Time 6.1367 -0.564
   Residual 4.9357       

Fixed effects: HR ~ Time + baseHR 
              Value Std.Error DF t-value p-value 
(Intercept)  36.931    9.6028 95  3.8459   2e-04
       Time  -7.027    1.8023 95 -3.8991   2e-04
     baseHR   0.551    0.1140 22  4.8325   1e-04
 Correlation: 
       (Intr)   Time 
  Time -0.100       
baseHR -0.984  0.000

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3    Max 
 -2.1929 -0.50894 -0.038365 0.42728 2.4281

Number of Observations: 120
Number of Groups: 24 
S> plot( augPred( fm4HR, length = 2 ) )  # compare predictions and data
