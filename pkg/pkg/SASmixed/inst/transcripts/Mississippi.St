S> ### $Id: Mississippi.St,v 1.1 1999/10/13 00:50:09 saikat Exp $
S> ### Analysis of the Mississippi nitrogren concentrations given as data set
S> ### 4.2 in "SAS System for Mixed Models"
S> options(contrasts=c(factor="contr.SAS", ordered="contr.poly"))
S> formula( Mississippi )
y ~ 1 | influent
S> plot( Mississippi )
S> fm1Miss <- lme( y ~ 1, data = Mississippi, random = ~ 1 | influent,
+                method = "ML")
S> summary( fm1Miss )         # compare with output 4.2, p. 143
Linear mixed-effects model fit by maximum likelihood
 Data: Mississippi 
     AIC    BIC  logLik 
  262.56 267.39 -128.28

Random effects:
 Formula:  ~ 1 | influent
        (Intercept) Residual 
StdDev:      7.1592   6.5343

Fixed effects: y ~ 1 
             Value Std.Error DF t-value p-value 
(Intercept) 21.217    3.1651 31  6.7035  <.0001

Standardized Within-Group Residuals:
     Min       Q1       Med     Q3    Max 
 -1.9344 -0.55704 -0.035454 0.8376 1.9323

Number of Observations: 37
Number of Groups: 6 
S> fm1RMiss <- update( fm1Miss, method = "REML" )
S> summary( fm1RMiss )        # compare with output 4.1, p. 142
Linear mixed-effects model fit by REML
 Data: Mississippi 
     AIC   BIC  logLik 
  258.35 263.1 -126.18

Random effects:
 Formula:  ~ 1 | influent
        (Intercept) Residual 
StdDev:      7.9576   6.5313

Fixed effects: y ~ 1 
             Value Std.Error DF t-value p-value 
(Intercept) 21.223     3.429 31  6.1892  <.0001

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3    Max 
 -1.9144 -0.53646 -0.032173 0.83714 1.9582

Number of Observations: 37
Number of Groups: 6 
S> random.effects( fm1Miss )    # BLUP's of random effects on p. 144
  (Intercept) 
3    -3.78627
5    -5.84352
2    -6.57723
1     0.30978
4     2.88267
6    13.01457
S> random.effects( fm1Miss , aug = TRUE )   # including covariates
  (Intercept)      y Type 
3    -3.78627 16.800    1
5    -5.84352 14.400    1
2    -6.57723 13.857    2
1     0.30978 21.556    2
4     2.88267 24.500    2
6    13.01457 36.400    3
S> plot( random.effects( fm1Miss , aug = TRUE ), outer = ~ Type )
Error in panel(x = c(-3.7862747633584, -5...: Can't pass data of
	mode "call" to old-S function
Dumped
S> random.effects( fm1RMiss )   # BLUP's of random effects on p. 142
  (Intercept) 
3    -3.89795
5    -6.01298
2    -6.71933
1     0.30929
4     2.94610
6    13.37487
S> intervals( fm1RMiss )        # interval estimates of variance components
Approximate 95% confidence intervals

 Fixed effects:
            lower   est.  upper 
(Intercept) 14.23 21.223 28.217

 Random Effects:
  Level: influent 
                 lower   est.  upper 
sd((Intercept)) 3.9513 7.9576 16.026

 Within-group standard error:
  lower   est.  upper 
 5.0895 6.5313 8.3815
S> c(2.9568, 7.9576, 21.416)^2 # compare to output 4.7, p. 148
[1]   8.7427  63.3234 458.6451
S> fm2RMiss <- lme( y ~ Type, data = Mississippi, random = ~ 1 | influent,
+       method = "REML" )
S> summary( fm2RMiss )         # compare to output 4.8 and 4.9, pp. 150-152
Linear mixed-effects model fit by REML
 Data: Mississippi 
     AIC    BIC  logLik 
  244.52 252.16 -117.26

Random effects:
 Formula:  ~ 1 | influent
        (Intercept) Residual 
StdDev:      3.8691   6.5202

Fixed effects: y ~ Type 
              Value Std.Error DF t-value p-value 
(Intercept)  36.400    4.8449 31  7.5131  <.0001
      Type1 -20.800    5.9338  3 -3.5054  0.0393
      Type2 -16.462    5.5168  3 -2.9840  0.0584
 Correlation: 
      (Intr)  Type1 
Type1 -0.816       
Type2 -0.878  0.717

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3    Max 
 -2.0877 -0.63918 -0.025703 0.70549 2.0191

Number of Observations: 37
Number of Groups: 6 
S> fm2Miss <- update( fm2RMiss, method = "ML" )  # get ML results too
S> anova( fm1Miss, fm2Miss )         # getting a p-value for the Type
        Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm1Miss     1  3 262.56 267.39 -128.28                       
fm2Miss     2  5 256.57 264.63 -123.29 1 vs 2  9.9834  0.0068
S>   ## Notice that the p-value is considerably smaller than for the F-test.
