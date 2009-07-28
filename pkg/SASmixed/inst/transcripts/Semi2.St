S> ### $Id: Semi2.St,v 1.1 1999/10/13 00:50:09 saikat Exp $
S> ### Analysis of the Oxide layer thickness data given as data set 4.4 in
S> ### "SAS System for Mixed Models"
S> options(contrasts = c(factor = "contr.SAS", ordered = "contr.poly"))
S> formula(Semi2)
Thickness ~ 1 | Lot/Wafer
S> plot(Semi2)
S> fm1Semi2 <- lme( Thickness ~ 1, data = Semi2,
+    random = ~ 1 | Lot/Wafer, method = "ML" )
S> summary( fm1Semi2 )
Linear mixed-effects model fit by maximum likelihood
 Data: Semi2 
     AIC    BIC  logLik 
  466.68 475.78 -229.34

Random effects:
 Formula:  ~ 1 | Lot
        (Intercept) 
StdDev:      10.583

 Formula:  ~ 1 | Wafer %in% Lot
        (Intercept) Residual 
StdDev:      5.9888   3.5453

Fixed effects: Thickness ~ 1 
             Value Std.Error DF t-value p-value 
(Intercept) 2000.2    3.9862 48  501.77  <.0001

Standardized Within-Group Residuals:
     Min       Q1     Med      Q3   Max 
 -1.8688 -0.49906 0.10125 0.54825 1.798

Number of Observations: 72
Number of Groups: 
 Lot Wafer %in% Lot 
   8             24
S> fm1RSemi2 <- update( fm1Semi2, method = "REML" )
S> summary( fm1RSemi2 )       # compare with output 4.13, p. 156
Linear mixed-effects model fit by REML
 Data: Semi2 
     AIC    BIC  logLik 
  462.02 471.07 -227.01

Random effects:
 Formula:  ~ 1 | Lot
        (Intercept) 
StdDev:      11.398

 Formula:  ~ 1 | Wafer %in% Lot
        (Intercept) Residual 
StdDev:      5.9888   3.5453

Fixed effects: Thickness ~ 1 
             Value Std.Error DF t-value p-value 
(Intercept) 2000.2    4.2317 48  472.66  <.0001

Standardized Within-Group Residuals:
     Min       Q1     Med      Q3    Max 
 -1.8746 -0.49908 0.10465 0.55102 1.7922

Number of Observations: 72
Number of Groups: 
 Lot Wafer %in% Lot 
   8             24
S> fm2Semi2 <- update( fm1Semi2, Thickness ~ Source )
S> anova( fm1Semi2, fm2Semi2 )
         Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm1Semi2     1  4 466.68 475.78 -229.34                       
fm2Semi2     2  5 466.86 478.25 -228.43 1 vs 2   1.813  0.1782
S> ## Again, the p-value is smaller than that for the F test.
S> fm2RSemi2 <- update( fm2Semi2, method = "REML" )
S> summary( fm2RSemi2 )       # compare with output 4.15, p. 159
Linear mixed-effects model fit by REML
 Data: Semi2 
     AIC    BIC  logLik 
  456.48 467.72 -223.24

Random effects:
 Formula:  ~ 1 | Lot
        (Intercept) 
StdDev:       10.95

 Formula:  ~ 1 | Wafer %in% Lot
        (Intercept) Residual 
StdDev:      5.9888   3.5453

Fixed effects: Thickness ~ Source 
              Value Std.Error DF t-value p-value 
(Intercept)  2005.2    5.7716 48  347.43  <.0001
     Source   -10.1    8.1622  6   -1.24  0.2629
 Correlation: 
       (Intr) 
Source -0.707

Standardized Within-Group Residuals:
     Min       Q1     Med      Q3    Max 
 -1.8864 -0.51397 0.11776 0.53466 1.7804

Number of Observations: 72
Number of Groups: 
 Lot Wafer %in% Lot 
   8             24
S> fm3Semi2 <- update( fm2Semi2, 
+    random = list(Lot = pdDiag( ~ Source - 1 ), Wafer = ~ 1 ) )
S> summary( fm3Semi2 )
Linear mixed-effects model fit by maximum likelihood
 Data: Semi2 
     AIC   BIC  logLik 
  465.24 478.9 -226.62

Random effects:
 Formula:  ~ Source - 1 | Lot
 Structure: Diagonal
        Source1 Source2 
StdDev:  3.0745  12.794

 Formula:  ~ 1 | Wafer %in% Lot
        (Intercept) Residual 
StdDev:      5.9891   3.5454

Fixed effects: Thickness ~ Source 
              Value Std.Error DF t-value p-value 
(Intercept)  2005.2    6.7474 48  297.18  <.0001
     Source   -10.1    7.1688  6   -1.41  0.2092
 Correlation: 
       (Intr) 
Source -0.941

Standardized Within-Group Residuals:
     Min       Q1      Med    Q3    Max 
 -1.8936 -0.55242 0.079282 0.575 1.7732

Number of Observations: 72
Number of Groups: 
 Lot Wafer %in% Lot 
   8             24
S> fm3RSemi2 <- update( fm3Semi2, method = "REML" )
S> summary( fm3RSemi2 )       # compare with output 4.17, p. 163
Linear mixed-effects model fit by REML
 Data: Semi2 
     AIC    BIC  logLik 
  455.76 469.25 -221.88

Random effects:
 Formula:  ~ Source - 1 | Lot
 Structure: Diagonal
        Source1 Source2 
StdDev:  4.1318  14.923

 Formula:  ~ 1 | Wafer %in% Lot
        (Intercept) Residual 
StdDev:      5.9889   3.5453

Fixed effects: Thickness ~ Source 
              Value Std.Error DF t-value p-value 
(Intercept)  2005.2    7.6821 48  261.02  <.0001
     Source   -10.1    8.1622  6   -1.24  0.2629
 Correlation: 
       (Intr) 
Source -0.941

Standardized Within-Group Residuals:
     Min       Q1      Med      Q3    Max 
 -1.8991 -0.53091 0.081684 0.55722 1.7677

Number of Observations: 72
Number of Groups: 
 Lot Wafer %in% Lot 
   8             24
S> anova( fm1Semi2, fm2Semi2, fm3Semi2 )
         Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm1Semi2     1  4 466.68 475.78 -229.34                       
fm2Semi2     2  5 466.86 478.25 -228.43 1 vs 2  1.8130  0.1782
fm3Semi2     3  6 465.24 478.90 -226.62 2 vs 3  3.6197  0.0571
