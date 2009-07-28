S> ### $Id: CultMiss.St,v 1.1 1999/10/13 00:50:09 saikat Exp $
S> ### A blocked split-plot with missing data (sec 2.7, pp. 68-75)
S> ## Remove all observations in block 1, cultivar 'A'.
S> options( contrasts = c(factor = "contr.SAS", ordered = "contr.poly") )
S> CultMiss <- Cultivation[ Cultivation$Block != 1 | Cultivation$Cult != 'a', ]
S> dim(CultMiss)
[1] 21  4
S> plot( CultMiss, inner = ~ Inoc )
S> fm1CultM <- lme( drywt ~ Cult * Inoc, CultMiss, list(Block = ~ 1, Cult = ~ 1),
+                 method = "ML")
S> summary( fm1CultM )
Linear mixed-effects model fit by maximum likelihood
 Data: CultMiss 
     AIC    BIC  logLik 
  76.647 86.047 -29.323

Random effects:
 Formula:  ~ 1 | Block
        (Intercept) 
StdDev:     0.86256

 Formula:  ~ 1 | Cult %in% Block
        (Intercept) Residual 
StdDev:     0.86879  0.66455

Fixed effects: drywt ~ Cult * Inoc 
              Value Std.Error DF t-value p-value 
(Intercept)  33.525    0.8241 10  40.680  <.0001
       Cult  -0.619    1.0123  2  -0.612  0.6031
    Inoccon  -5.500    0.5560 10  -9.892  <.0001
    Inocdea  -2.875    0.5560 10  -5.171  0.0004
CultInoccon   0.867    0.8493 10   1.020  0.3316
CultInocdea  -0.725    0.8493 10  -0.854  0.4133
 Correlation: 
            (Intr)   Cult Inoccn Inocde CltIncc 
       Cult -0.502                             
    Inoccon -0.337  0.275                      
    Inocdea -0.337  0.275  0.500               
CultInoccon  0.221 -0.420 -0.655 -0.327        
CultInocdea  0.221 -0.420 -0.327 -0.655  0.500 

Standardized Within-Group Residuals:
     Min       Q1     Med      Q3    Max 
 -1.4566 -0.55593 0.19646 0.37888 1.9535

Number of Observations: 21
Number of Groups: 
 Block Cult %in% Block 
     4               7
S> fm2CultM <- update( fm1CultM, drywt ~ Cult + Inoc )
S> fm3CultM <- update( fm1CultM, drywt ~ Inoc )
S> fm4CultM <- update( fm1CultM, drywt ~ 1 )
S> anova( fm1Cult, fm2Cult, fm3Cult, fm4Cult )
        Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm1Cult     1  9  86.76  97.36 -34.378                       
fm2Cult     2  7  85.88  94.13 -35.940 1 vs 2   3.125  0.2096
fm3Cult     3  6  84.79  91.85 -36.393 2 vs 3   0.905  0.3415
fm4Cult     4  4 122.48 127.19 -57.238 3 vs 4  41.691  <.0001
S> # Essentially the same conclusions as for the balanced data
S> fm3RCultM <- update( fm3CultM, method = "REML" )
S> summary( fm3RCultM )
Linear mixed-effects model fit by REML
 Data: CultMiss 
     AIC    BIC  logLik 
  74.474 79.816 -31.237

Random effects:
 Formula:  ~ 1 | Block
        (Intercept) 
StdDev:      1.1726

 Formula:  ~ 1 | Cult %in% Block
        (Intercept) Residual 
StdDev:     0.83473  0.83465

Fixed effects: drywt ~ Inoc 
              Value Std.Error DF t-value p-value 
(Intercept)  33.305   0.74204 12  44.883  <.0001
    Inoccon  -5.129   0.44614 12 -11.495  <.0001
    Inocdea  -3.186   0.44614 12  -7.141  <.0001
 Correlation: 
        (Intr) Inoccn 
Inoccon -0.301       
Inocdea -0.301  0.500

Standardized Within-Group Residuals:
     Min      Q1      Med      Q3    Max 
 -1.7182 -0.2666 -0.14867 0.36054 1.7434

Number of Observations: 21
Number of Groups: 
 Block Cult %in% Block 
     4               7
