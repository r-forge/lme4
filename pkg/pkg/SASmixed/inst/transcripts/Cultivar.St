S> ### $Id: Cultivar.St,v 1.1 1999/10/13 00:50:09 saikat Exp $
S> ### Split-plot experiment with whole plots in randomized blocks (sec. 2.5)
S> options( contrasts = c(factor = "contr.SAS", ordered = "contr.poly") )
S> names( Cultivation )
[1] "Block" "Cult"  "Inoc"  "drywt"
S> formula( Cultivation )
drywt ~ 1 | Block/Cult
S> plot( Cultivation, inner = ~ Inoc )
S> fm1Cult <- lme( drywt ~ Inoc * Cult, 
+                  data = Cultivation, method = "ML",
+                  random = list( Block = ~ 1, Cult = ~ 1 ) )
S> summary( fm1Cult )
Linear mixed-effects model fit by maximum likelihood
 Data: Cultivation 
     AIC    BIC  logLik 
  86.755 97.358 -34.378

Random effects:
 Formula:  ~ 1 | Block
        (Intercept) 
StdDev:     0.81248

 Formula:  ~ 1 | Cult %in% Block
        (Intercept) Residual 
StdDev:     0.78329  0.72737

Fixed effects: drywt ~ Inoc * Cult 
              Value Std.Error DF t-value p-value 
(Intercept)  33.525   0.77518 12  43.248  <.0001
    Inoccon  -5.500   0.59390 12  -9.261  <.0001
    Inocdea  -2.875   0.59390 12  -4.841  0.0004
       Cult  -0.375   0.87278  3  -0.430  0.6964
InocconCult   0.250   0.83990 12   0.298  0.7711
InocdeaCult  -1.025   0.83990 12  -1.220  0.2458
 Correlation: 
            (Intr) Inoccn Inocde   Cult InccnC 
    Inoccon -0.383                            
    Inocdea -0.383  0.500                     
       Cult -0.563  0.340  0.340              
InocconCult  0.271 -0.707 -0.354 -0.481       
InocdeaCult  0.271 -0.354 -0.707 -0.481  0.500

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3   Max 
 -1.4244 -0.54855 -0.018575 0.73293 1.874

Number of Observations: 24
Number of Groups: 
 Block Cult %in% Block 
     4               8
S> logLik( update( fm1Cult, method = "REML" ) )  # check the same model is being fit
[1] -32.531
S> fm2Cult <- update( fm1Cult, drywt ~ Inoc + Cult )
S> fm3Cult <- update( fm1Cult, drywt ~ Inoc )
S> fm4Cult <- update( fm1Cult, drywt ~ 1 )
S> anova( fm1Cult, fm2Cult, fm3Cult, fm4Cult )
        Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm1Cult     1  9  86.76  97.36 -34.378                       
fm2Cult     2  7  85.88  94.13 -35.940 1 vs 2   3.125  0.2096
fm3Cult     3  6  84.79  91.85 -36.393 2 vs 3   0.905  0.3415
fm4Cult     4  4 122.48 127.19 -57.238 3 vs 4  41.691  <.0001
S> fm5Cult <- update( fm1Cult, drywt ~ Cult )
S> anova( fm1Cult, fm2Cult, fm5Cult, fm4Cult )
        Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm1Cult     1  9  86.76  97.36 -34.378                       
fm2Cult     2  7  85.88  94.13 -35.940 1 vs 2   3.125  0.2096
fm5Cult     3  5 124.12 130.02 -57.062 2 vs 3  42.245  <.0001
fm4Cult     4  4 122.48 127.19 -57.238 3 vs 4   0.351  0.5534
S> ### AIC, BIC, and likelihood ratio tests all prefer fm3Cult
S> summary( update( fm3Cult, method = "REML" ) )  # REML estimates
Linear mixed-effects model fit by REML
 Data: Cultivation 
    AIC    BIC logLik 
  84.28 90.547 -36.14

Random effects:
 Formula:  ~ 1 | Block
        (Intercept) 
StdDev:     0.95452

 Formula:  ~ 1 | Cult %in% Block
        (Intercept) Residual 
StdDev:     0.86372  0.85731

Fixed effects: drywt ~ Inoc 
              Value Std.Error DF t-value p-value 
(Intercept)  33.337   0.64257 14  51.881  <.0001
    Inoccon  -5.375   0.42865 14 -12.539  <.0001
    Inocdea  -3.388   0.42865 14  -7.903  <.0001
 Correlation: 
        (Intr) Inoccn 
Inoccon -0.334       
Inocdea -0.334  0.500

Standardized Within-Group Residuals:
     Min       Q1      Med      Q3    Max 
 -1.5717 -0.41039 -0.08165 0.53191 1.5561

Number of Observations: 24
Number of Groups: 
 Block Cult %in% Block 
     4               8
