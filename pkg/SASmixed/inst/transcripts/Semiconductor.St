S> ### Analysis of the Semiconductor split-plot experiment, section 2.3 of
S> ### "SAS System for Mixed Models"
S> options(contrasts = c(factor = "contr.SAS", ordered = "contr.poly"))
S> formula( Semiconductor )
resistance ~ 1 | Grp
S> plot( Semiconductor )
S> names( Semiconductor )
[1] "resistance" "ET"         "Wafer"      "position"   "Grp"       
S> plot( Semiconductor, inner = ~ position )
S> fm1Semi <- lme( resistance ~ ET * position,
+                  data = Semiconductor,
+                  random = ~ 1 | Grp, method = "ML" )
S> summary( fm1Semi )
Linear mixed-effects model fit by maximum likelihood
 Data: Semiconductor 
     AIC    BIC  logLik 
  66.147 99.828 -15.073

Random effects:
 Formula:  ~ 1 | Grp
        (Intercept) Residual 
StdDev:     0.26557  0.27221

Fixed effects: resistance ~ ET * position 
               Value Std.Error DF t-value p-value 
 (Intercept)  6.5400   0.26891 24  24.320  <.0001
         ET1 -0.6533   0.38030  8  -1.718  0.1241
         ET2 -0.6233   0.38030  8  -1.639  0.1398
         ET3 -0.4467   0.38030  8  -1.175  0.2740
   position1 -0.2000   0.27221 24  -0.735  0.4696
   position2  0.0133   0.27221 24   0.049  0.9613
   position3 -0.6433   0.27221 24  -2.363  0.0266
ET1position1 -0.0733   0.38497 24  -0.190  0.8505
ET2position1  0.2767   0.38497 24   0.719  0.4793
ET3position1  0.2433   0.38497 24   0.632  0.5333
ET1position2 -0.4500   0.38497 24  -1.169  0.2539
ET2position2  0.2567   0.38497 24   0.667  0.5113
ET3position2  0.2400   0.38497 24   0.623  0.5389
ET1position3  0.3100   0.38497 24   0.805  0.4286
ET2position3  0.4933   0.38497 24   1.281  0.2123
ET3position3  0.3233   0.38497 24   0.840  0.4093
 Correlation: 
             (Intr)    ET1    ET2    ET3 postn1 postn2 postn3 
         ET1 -0.707                                          
         ET2 -0.707  0.500                                   
         ET3 -0.707  0.500  0.500                            
   position1 -0.506  0.358  0.358  0.358                     
   position2 -0.506  0.358  0.358  0.358  0.500              
   position3 -0.506  0.358  0.358  0.358  0.500  0.500       
ET1position1  0.358 -0.506 -0.253 -0.253 -0.707 -0.354 -0.354
ET2position1  0.358 -0.253 -0.506 -0.253 -0.707 -0.354 -0.354
ET3position1  0.358 -0.253 -0.253 -0.506 -0.707 -0.354 -0.354
ET1position2  0.358 -0.506 -0.253 -0.253 -0.354 -0.707 -0.354
ET2position2  0.358 -0.253 -0.506 -0.253 -0.354 -0.707 -0.354
ET3position2  0.358 -0.253 -0.253 -0.506 -0.354 -0.707 -0.354
ET1position3  0.358 -0.506 -0.253 -0.253 -0.354 -0.354 -0.707
ET2position3  0.358 -0.253 -0.506 -0.253 -0.354 -0.354 -0.707
ET3position3  0.358 -0.253 -0.253 -0.506 -0.354 -0.354 -0.707
             ET1ps1 ET2ps1 ET3ps1 ET1ps2 ET2ps2 ET3ps2 ET1ps3 
         ET1                                                 
         ET2                                                 
         ET3                                                 
   position1                                                 
   position2                                                 
   position3                                                 
ET1position1                                                 
ET2position1  0.500                                          
ET3position1  0.500  0.500                                   
ET1position2  0.500  0.250  0.250                            
ET2position2  0.250  0.500  0.250  0.500                     
ET3position2  0.250  0.250  0.500  0.500  0.500              
ET1position3  0.500  0.250  0.250  0.500  0.250  0.250       
ET2position3  0.250  0.500  0.250  0.250  0.500  0.250  0.500
ET3position3  0.250  0.250  0.500  0.250  0.250  0.500  0.500
             ET2ps3 
         ET1       
         ET2       
         ET3       
   position1       
   position2       
   position3       
ET1position1       
ET2position1       
ET3position1       
ET1position2       
ET2position2       
ET3position2       
ET1position3       
ET2position3       
ET3position3  0.500

Standardized Within-Group Residuals:
     Min      Q1    Med      Q3    Max 
 -2.3406 -0.5624 0.0126 0.57402 1.6062

Number of Observations: 48
Number of Groups: 12 
S> ## check significance of the interaction
S> fm2Semi <- update( fm1Semi, resistance ~ ET + position )
S> fm3Semi <- update( fm1Semi, resistance ~ ET )
S> fm4Semi <- update( fm1Semi, resistance ~ position )
S> fm5Semi <- update( fm1Semi, resistance ~ 1 )
S> anova( fm5Semi, fm4Semi, fm2Semi, fm1Semi )
        Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm5Semi     1  3 62.373 67.987 -28.187                       
fm4Semi     2  6 58.251 69.479 -23.126 1 vs 2  10.122  0.0176
fm2Semi     3  9 57.687 74.528 -19.844 2 vs 3   6.564  0.0872
fm1Semi     4 18 66.147 99.828 -15.073 3 vs 4   9.540  0.3889
S> anova( fm5Semi, fm3Semi, fm2Semi, fm1Semi )
        Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm5Semi     1  3 62.373 67.987 -28.187                       
fm3Semi     2  6 61.809 73.036 -24.905 1 vs 2   6.564  0.0872
fm2Semi     3  9 57.687 74.528 -19.844 2 vs 3  10.122  0.0176
fm1Semi     4 18 66.147 99.828 -15.073 3 vs 4   9.540  0.3889
S> ## AIC favors resistance ~ ET + position
S> ## BIC favors resistance ~ 1
S> ## Likelihood ratio seems to favor resistance ~ position
