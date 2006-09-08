S> ### $Id: Bond.St,v 1.1 1999/10/13 00:50:09 saikat Exp $
S> ## Set options 
S> options( prompt = "S> ", digits = 5, width = 65, 
+   contrasts = c(factor = "contr.SAS", ordered = "contr.poly") )
S> plot(Bond)                       # dotplot by Ingot
S> plot(Bond, inner = ~ Metal)      # different symbols for different Metals
S> formula(Bond)                    # check the formula
pressure ~ 1 | Ingot
S> fm1Bond <- lme( pressure ~ Metal, data = Bond, random = ~ 1 | Ingot,
+                method = "ML")
S> summary( fm1Bond )
Linear mixed-effects model fit by maximum likelihood
 Data: Bond 
     AIC    BIC  logLik 
  125.71 130.93 -57.854

Random effects:
 Formula:  ~ 1 | Ingot
        (Intercept) Residual 
StdDev:      3.1325   2.9816

Fixed effects: pressure ~ Metal 
              Value Std.Error DF t-value p-value 
(Intercept)  71.100    1.7655 12  40.271  <.0001
     Metalc  -0.914    1.7214 12  -0.531  0.6050
     Metali   4.800    1.7214 12   2.788  0.0164
 Correlation: 
       (Intr) Metalc 
Metalc -0.488       
Metali -0.488  0.500

Standardized Within-Group Residuals:
     Min       Q1     Med      Q3    Max 
 -1.4551 -0.81901 0.08048 0.52228 1.9611

Number of Observations: 21
Number of Groups: 7 
S> ## default criterion in lme is maximum likelihood (ML).
S> ## Re-fit to get REML results
S> fm1RBond <- update( fm1Bond, method = "REML" )
S> summary( fm1RBond )              # compare with output 1.1 on p. 6
Linear mixed-effects model fit by REML
 Data: Bond 
     AIC    BIC  logLik 
  117.79 122.24 -53.895

Random effects:
 Formula:  ~ 1 | Ingot
        (Intercept) Residual 
StdDev:      3.3835   3.2205

Fixed effects: pressure ~ Metal 
              Value Std.Error DF t-value p-value 
(Intercept)  71.100    1.7655 12  40.271  <.0001
     Metalc  -0.914    1.7214 12  -0.531  0.6050
     Metali   4.800    1.7214 12   2.788  0.0164
 Correlation: 
       (Intr) Metalc 
Metalc -0.488       
Metali -0.488  0.500

Standardized Within-Group Residuals:
     Min       Q1     Med      Q3    Max 
 -1.3471 -0.75825 0.07451 0.48354 1.8157

Number of Observations: 21
Number of Groups: 7 
S> logLik( fm1RBond )                # log-restricted-likelihood
[1] -53.895
S> c(3.3835,3.2205)^2                # variance estimates 
[1] 11.448 10.372
S> ##  To test the need for the Metal term in the fixed effects, 
S> ##  re-fit and use anova.  You must use ML to do this.  RML results are
S> ##  not comparable if the fixed-effects specification changes.
S> fm2Bond <- update( fm1Bond, pressure ~ 1 )
S> anova( fm1Bond, fm2Bond )
        Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm1Bond     1  5 125.71 130.93 -57.854                       
fm2Bond     2  3 131.82 134.96 -62.912 1 vs 2  10.116  0.0064
