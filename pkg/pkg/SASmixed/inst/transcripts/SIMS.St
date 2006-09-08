S> ### $Id: SIMS.St,v 1.1 1999/10/13 00:50:09 saikat Exp $
S> ### Analysis of the Second International Mathematics Study (SIMS)
S> ### described in section 7.2 of "SAS System for Mixed Models"
S> plot(SIMS, layout = c(10, 10, 2) )  # high level of variability within classes
S> ## to record the execution time for comparison the fit was done as
S> ## unix.time(assign("fm1RSIMS", lme(Gain ~ Pretot, SIMS, ~ Pretot | Class, REML = TRUE, control = list(msVerbose = TRUE))))
S> ## Since unix.time is not available on Windows, we give the equivalent
S> ## statement here
S> fm1RSIMS <- lme(Gain ~ Pretot, data = SIMS, 
+     random = ~ Pretot | Class, control = list(msVerbose = TRUE))
Iteration:  0 ,  1  function calls, F=  21105 
Parameters:
[1] 0.205678 0.066285 4.012905
Iteration:  1 ,  2  function calls, F=  21105 
Parameters:
[1] 0.213001 0.064138 4.126911
Iteration:  2 ,  3  function calls, F=  21105 
Parameters:
[1] 0.215042 0.063723 4.156579
S> ## Timing done on an UltraSPARC 2 running S-PLUS 3.4 for Solaris.
S> ## Your mileage may vary.
S> summary(fm1RSIMS)              # compare to output 7.4, p. 262
Linear mixed-effects model fit by REML
 Data: SIMS 
    AIC   BIC logLik 
  22393 22430 -11190

Random effects:
 Formula:  ~ Pretot | Class
 Structure: General positive-definite
              StdDev   Corr 
(Intercept) 3.806502 (Inter
     Pretot 0.095932 -0.641
   Residual 4.715480       

Fixed effects: Gain ~ Pretot 
              Value Std.Error   DF t-value p-value 
(Intercept)  7.0596    0.3659 3500  19.294  <.0001
     Pretot -0.1860    0.0161 3500 -11.556  <.0001
 Correlation: 
       (Intr) 
Pretot -0.76 

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3    Max 
 -5.3776 -0.59554 0.0098766 0.64099 3.7183

Number of Observations: 3691
Number of Groups: 190 
S> intervals( fm1RSIMS )
Approximate 95% confidence intervals

 Fixed effects:
              lower     est.    upper 
(Intercept)  6.3422  7.05961  7.77700
     Pretot -0.2176 -0.18603 -0.15447

 Random Effects:
  Level: Class 
                            lower      est.    upper 
        sd((Intercept))  3.172870  3.806502  4.56667
             sd(Pretot)  0.056691  0.095932  0.16233
cor((Intercept),Pretot) -0.834252 -0.640782 -0.30683

 Within-group standard error:
  lower   est.  upper 
 4.6022 4.7155 4.8316
