S> ### $Id: Weights.St,v 1.1 1999/10/13 00:50:09 saikat Exp $
S> ### $Id: Weights.St,v 1.1 1999/10/13 00:50:09 saikat Exp $
S> ### Analysis of the weight-lifting program data given as data set 3.2(a)
S> ### in "SAS System for Mixed Models"
S> options(contrasts = c(factor = "contr.SAS", ordered = "contr.poly"))
S> plot(Weights, layout = c(21,3), skip = rep(c(F,T,F,T,F),c(20,1,16,5,21)))
S> fm1Weight <- 
+    lme( strength ~ Program * Time, data = Weights, random = ~ 1 | Subj,
+        method = "ML" )
S> summary( fm1Weight )
Linear mixed-effects model fit by maximum likelihood
 Data: Weights 
     AIC    BIC  logLik 
  1441.9 1473.8 -712.93

Random effects:
 Formula:  ~ 1 | Subj
        (Intercept) Residual 
StdDev:      3.0152    1.085

Fixed effects: strength ~ Program * Time 
                  Value Std.Error  DF t-value p-value 
    (Intercept)  81.105    0.6874 339  117.98  <.0001
    ProgramCONT  -1.115    0.9843  54   -1.13  0.2622
      ProgramRI  -1.045    1.0454  54   -1.00  0.3219
           Time   0.160    0.0225 339    7.09  <.0001
ProgramCONTTime  -0.184    0.0323 339   -5.70  <.0001
  ProgramRITime  -0.055    0.0343 339   -1.60  0.1098
 Correlation: 
                (Intr) PrCONT PrgrRI   Time PCONTT 
    ProgramCONT -0.698                            
      ProgramRI -0.658  0.459                     
           Time -0.230  0.160  0.151              
ProgramCONTTime  0.160 -0.230 -0.105 -0.698       
  ProgramRITime  0.151 -0.105 -0.230 -0.658  0.459

Standardized Within-Group Residuals:
     Min       Q1      Med      Q3    Max 
 -3.1292 -0.61237 0.037615 0.63314 3.3469

Number of Observations: 399
Number of Groups: 57 
S> summary( update( fm1Weight, method = "REML" ) )  # compare with output 3.1, p. 91
Linear mixed-effects model fit by REML
 Data: Weights 
     AIC    BIC  logLik 
  1455.4 1487.2 -719.68

Random effects:
 Formula:  ~ 1 | Subj
        (Intercept) Residual 
StdDev:      3.0991   1.0897

Fixed effects: strength ~ Program * Time 
                  Value Std.Error  DF t-value p-value 
    (Intercept)  81.105    0.7001 339  115.84  <.0001
    ProgramCONT  -1.115    1.0024  54   -1.11  0.2708
      ProgramRI  -1.045    1.0647  54   -0.98  0.3306
           Time   0.160    0.0225 339    7.11  <.0001
ProgramCONTTime  -0.184    0.0322 339   -5.72  <.0001
  ProgramRITime  -0.055    0.0342 339   -1.61  0.1087
 Correlation: 
                (Intr) PrCONT PrgrRI   Time PCONTT 
    ProgramCONT -0.698                            
      ProgramRI -0.658  0.459                     
           Time -0.225  0.157  0.148              
ProgramCONTTime  0.157 -0.225 -0.103 -0.698       
  ProgramRITime  0.148 -0.103 -0.225 -0.658  0.459

Standardized Within-Group Residuals:
     Min       Q1      Med      Q3    Max 
 -3.1167 -0.60963 0.039635 0.63075 3.3352

Number of Observations: 399
Number of Groups: 57 
S> c( 3.0991, 1.0897 )^2
[1] 9.6044 1.1874
S> fm2Weight <- update( fm1Weight, random = ~ Time | Subj )
S> anova( fm1Weight, fm2Weight )
          Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm1Weight     1  8 1441.9 1473.8 -712.93                       
fm2Weight     2 10 1333.4 1373.2 -656.68 1 vs 2  112.52  <.0001
S> plot(augPred( fm2Weight ), layout = c(21,3),
+       skip = rep(c(F,T,F,T,F),c(20,1,16,5,21)))
S> summary( fm2Weight )
Linear mixed-effects model fit by maximum likelihood
 Data: Weights 
     AIC    BIC  logLik 
  1333.4 1373.2 -656.68

Random effects:
 Formula:  ~ Time | Subj
 Structure: General positive-definite
             StdDev   Corr 
(Intercept) 2.92292 (Inter
       Time 0.17074 -0.115
   Residual 0.79559       

Fixed effects: strength ~ Program * Time 
                  Value Std.Error  DF t-value p-value 
    (Intercept)  81.105   0.65635 339  123.57  <.0001
    ProgramCONT  -1.115   0.93976  54   -1.19  0.2405
      ProgramRI  -1.045   0.99811  54   -1.05  0.2997
           Time   0.160   0.04102 339    3.90  0.0001
ProgramCONTTime  -0.184   0.05873 339   -3.13  0.0019
  ProgramRITime  -0.055   0.06238 339   -0.88  0.3790
 Correlation: 
                (Intr) PrCONT PrgrRI   Time PCONTT 
    ProgramCONT -0.698                            
      ProgramRI -0.658  0.459                     
           Time -0.174  0.122  0.114              
ProgramCONTTime  0.122 -0.174 -0.080 -0.698       
  ProgramRITime  0.114 -0.080 -0.174 -0.658  0.459

Standardized Within-Group Residuals:
     Min      Q1       Med      Q3    Max 
 -2.7527 -0.6483 -0.010436 0.57345 3.4828

Number of Observations: 399
Number of Groups: 57 
S> fm3Weight <- update( fm2Weight, correlation = corAR1())
S> anova( fm2Weight, fm3Weight )
          Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm2Weight     1 10 1333.4 1373.2 -656.68                       
fm3Weight     2 11 1284.2 1328.1 -631.12 1 vs 2  51.112  <.0001
S> fm4Weight <- update( fm3Weight, strength ~ Program * (Time + I(Time^2)),
+                     random = ~Time|Subj)
S> anova( fm1Weight, fm2Weight, fm3Weight, fm4Weight )
          Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm1Weight     1  8 1441.9 1473.8 -712.93                       
fm2Weight     2 10 1333.4 1373.2 -656.68 1 vs 2  112.52  <.0001
fm3Weight     3 11 1284.2 1328.1 -631.12 2 vs 3   51.11  <.0001
fm4Weight     4 14 1281.1 1336.9 -626.53 3 vs 4    9.18   0.027
S> summary( fm4Weight )
Linear mixed-effects model fit by maximum likelihood
 Data: Weights 
     AIC    BIC  logLik 
  1281.1 1336.9 -626.53

Random effects:
 Formula:  ~ Time | Subj
 Structure: General positive-definite
              StdDev   Corr 
(Intercept) 2.441731 (Inter
       Time 0.081467 0.525 
   Residual 1.516788       

Correlation Structure: AR(1)
 Parameter estimate(s):
     Phi 
 0.78658
Fixed effects: strength ~ Program + Time + I(Time^2) + Program:Time + Program:I( 	Time^2) 
                       Value Std.Error  DF t-value p-value 
         (Intercept)  80.787    0.6601 336  122.38  <.0001
         ProgramCONT  -1.109    0.9451  54   -1.17  0.2456
           ProgramRI  -1.455    1.0038  54   -1.45  0.1530
                Time   0.277    0.0991 336    2.80  0.0055
           I(Time^2)  -0.008    0.0065 336   -1.19  0.2366
     ProgramCONTTime  -0.188    0.1419 336   -1.32  0.1864
       ProgramRITime   0.111    0.1507 336    0.74  0.4606
ProgramCONTI(Time^2)   0.000    0.0093 336    0.04  0.9718
  ProgramRII(Time^2)  -0.012    0.0099 336   -1.17  0.2413
 Correlation: 
                     (Intr) PrCONT PrgrRI   Time I(T^2) PCONTT 
         ProgramCONT -0.698                                   
           ProgramRI -0.658  0.459                            
                Time -0.255  0.178  0.168                     
           I(Time^2)  0.224 -0.157 -0.148 -0.920              
     ProgramCONTTime  0.178 -0.255 -0.117 -0.698  0.643       
       ProgramRITime  0.168 -0.117 -0.255 -0.658  0.605  0.459
ProgramCONTI(Time^2) -0.157  0.224  0.103  0.643 -0.698 -0.920
  ProgramRII(Time^2) -0.148  0.103  0.224  0.605 -0.658 -0.423
                     PrgRIT PCONTI(T^2) 
         ProgramCONT                   
           ProgramRI                   
                Time                   
           I(Time^2)                   
     ProgramCONTTime                   
       ProgramRITime                   
ProgramCONTI(Time^2) -0.423            
  ProgramRII(Time^2) -0.920  0.459     

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3    Max 
 -1.7367 -0.49145 -0.045068 0.47571 2.3992

Number of Observations: 399
Number of Groups: 57 
S> intervals( fm4Weight )
Approximate 95% confidence intervals

 Fixed effects:
                          lower        est.      upper 
         (Intercept)  79.502852  80.7865923 82.0703328
         ProgramCONT  -2.982768  -1.1093822  0.7640032
           ProgramRI  -3.444677  -1.4549607  0.5347552
                Time   0.084490   0.2772186  0.4699476
           I(Time^2)  -0.020387  -0.0077219  0.0049433
     ProgramCONTTime  -0.463819  -0.1878731  0.0880728
       ProgramRITime  -0.181751   0.1113306  0.4044118
ProgramCONTI(Time^2)  -0.017804   0.0003300  0.0184639
  ProgramRII(Time^2)  -0.030885  -0.0116254  0.0076346

 Random Effects:
  Level: Subj 
                           lower     est.  upper 
      sd((Intercept))  1.2718625 2.441731 4.6877
             sd(Time)  0.0048644 0.081467 1.3644
cor((Intercept),Time) -0.9989554 0.524765 0.9999

 Correlation structure:
       lower    est.  upper 
Phi -0.32307 0.78658 0.9855

 Within-group standard error:
   lower   est.  upper 
 0.40432 1.5168 5.6901
