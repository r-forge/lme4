S> ### $Id: Multilocation.St,v 1.1 1999/10/13 00:50:09 saikat Exp $
S> ### Analysis of the Multilocation data with fixed effects for the locations
S> options( contrasts = c(factor = "contr.SAS", ordered = "contr.poly") )
S> formula( Multilocation )
Adj ~ 1 | Location/Block
S> names( Multilocation )
[1] "obs"      "Location" "Block"    "Trt"      "Adj"     
[6] "Fe"       "Grp"     
S> ### Create a Block %in% Location factor
S> Multilocation$Grp <-
+   getGroups( Multilocation, form = ~ Location/Block, level = 2 )
S> fm1Mult <- lme( Adj ~ Location * Trt, data = Multilocation, ~ 1 | Grp,
+                method = "ML")
S> summary( fm1Mult )
Linear mixed-effects model fit by maximum likelihood
 Data: Multilocation 
      AIC    BIC logLik 
  -11.146 90.775 43.573

Random effects:
 Formula:  ~ 1 | Grp
        (Intercept) Residual 
StdDev:    0.061206  0.15183

Fixed effects: Adj ~ Location * Trt 
                   Value Std.Error DF t-value p-value 
     (Intercept)  2.8657   0.03859 54  74.269  <.0001
      Location.L  0.8136   0.11576 18   7.029  <.0001
      Location.Q  0.2440   0.11576 18   2.108  0.0493
      Location.C  0.0990   0.11576 18   0.855  0.4036
    Location ^ 4 -0.2996   0.11576 18  -2.588  0.0186
    Location ^ 5  0.1272   0.11576 18   1.099  0.2863
    Location ^ 6 -0.0331   0.11576 18  -0.286  0.7785
    Location ^ 7 -0.4080   0.11576 18  -3.525  0.0024
    Location ^ 8  0.0851   0.11576 18   0.735  0.4719
            Trt1  0.0583   0.05061 54   1.153  0.2541
            Trt2 -0.1880   0.05061 54  -3.715  0.0005
            Trt3  0.0838   0.05061 54   1.656  0.1036
  Location.LTrt1  0.2002   0.15183 54   1.319  0.1928
  Location.QTrt1 -0.3034   0.15183 54  -1.999  0.0507
  Location.CTrt1 -0.0409   0.15183 54  -0.269  0.7886
Location ^ 4Trt1  0.3640   0.15183 54   2.397  0.0200
Location ^ 5Trt1 -0.1480   0.15183 54  -0.975  0.3340
Location ^ 6Trt1  0.0246   0.15183 54   0.162  0.8718
Location ^ 7Trt1  0.1470   0.15183 54   0.968  0.3374
Location ^ 8Trt1  0.1968   0.15183 54   1.296  0.2004
  Location.LTrt2  0.0812   0.15183 54   0.535  0.5949
  Location.QTrt2 -0.3237   0.15183 54  -2.132  0.0376
  Location.CTrt2 -0.1377   0.15183 54  -0.907  0.3684
Location ^ 4Trt2  0.0713   0.15183 54   0.470  0.6403
Location ^ 5Trt2  0.0396   0.15183 54   0.261  0.7953
Location ^ 6Trt2  0.1124   0.15183 54   0.740  0.4624
Location ^ 7Trt2  0.2366   0.15183 54   1.558  0.1250
Location ^ 8Trt2 -0.0375   0.15183 54  -0.247  0.8060
  Location.LTrt3  0.1344   0.15183 54   0.885  0.3799
  Location.QTrt3 -0.3384   0.15183 54  -2.229  0.0300
  Location.CTrt3  0.0969   0.15183 54   0.638  0.5259
Location ^ 4Trt3  0.1093   0.15183 54   0.720  0.4746
Location ^ 5Trt3 -0.1372   0.15183 54  -0.904  0.3702
Location ^ 6Trt3  0.1836   0.15183 54   1.209  0.2319
Location ^ 7Trt3  0.2799   0.15183 54   1.844  0.0707
Location ^ 8Trt3 -0.0057   0.15183 54  -0.038  0.9700
 Correlation: 
                 (Intr) Lctn.L Lctn.Q Lctn.C Lctn^4 Lctn^5 
      Location.L  0.000                                   
      Location.Q  0.000  0.000                            
      Location.C  0.000  0.000  0.000                     
    Location ^ 4  0.000  0.000  0.000  0.000              
    Location ^ 5  0.000  0.000  0.000  0.000  0.000       
    Location ^ 6  0.000  0.000  0.000  0.000  0.000  0.000
    Location ^ 7  0.000  0.000  0.000  0.000  0.000  0.000
    Location ^ 8  0.000  0.000  0.000  0.000  0.000  0.000
            Trt1 -0.656  0.000  0.000  0.000  0.000  0.000
            Trt2 -0.656  0.000  0.000  0.000  0.000  0.000
            Trt3 -0.656  0.000  0.000  0.000  0.000  0.000
  Location.LTrt1  0.000 -0.656  0.000  0.000  0.000  0.000
  Location.QTrt1  0.000  0.000 -0.656  0.000  0.000  0.000
  Location.CTrt1  0.000  0.000  0.000 -0.656  0.000  0.000
Location ^ 4Trt1  0.000  0.000  0.000  0.000 -0.656  0.000
Location ^ 5Trt1  0.000  0.000  0.000  0.000  0.000 -0.656
Location ^ 6Trt1  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 7Trt1  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 8Trt1  0.000  0.000  0.000  0.000  0.000  0.000
  Location.LTrt2  0.000 -0.656  0.000  0.000  0.000  0.000
  Location.QTrt2  0.000  0.000 -0.656  0.000  0.000  0.000
  Location.CTrt2  0.000  0.000  0.000 -0.656  0.000  0.000
Location ^ 4Trt2  0.000  0.000  0.000  0.000 -0.656  0.000
Location ^ 5Trt2  0.000  0.000  0.000  0.000  0.000 -0.656
Location ^ 6Trt2  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 7Trt2  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 8Trt2  0.000  0.000  0.000  0.000  0.000  0.000
  Location.LTrt3  0.000 -0.656  0.000  0.000  0.000  0.000
  Location.QTrt3  0.000  0.000 -0.656  0.000  0.000  0.000
  Location.CTrt3  0.000  0.000  0.000 -0.656  0.000  0.000
Location ^ 4Trt3  0.000  0.000  0.000  0.000 -0.656  0.000
Location ^ 5Trt3  0.000  0.000  0.000  0.000  0.000 -0.656
Location ^ 6Trt3  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 7Trt3  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 8Trt3  0.000  0.000  0.000  0.000  0.000  0.000
                 Lctn^6 Lctn^7 Lctn^8   Trt1   Trt2   Trt3 
      Location.L                                          
      Location.Q                                          
      Location.C                                          
    Location ^ 4                                          
    Location ^ 5                                          
    Location ^ 6                                          
    Location ^ 7  0.000                                   
    Location ^ 8  0.000  0.000                            
            Trt1  0.000  0.000  0.000                     
            Trt2  0.000  0.000  0.000  0.500              
            Trt3  0.000  0.000  0.000  0.500  0.500       
  Location.LTrt1  0.000  0.000  0.000  0.000  0.000  0.000
  Location.QTrt1  0.000  0.000  0.000  0.000  0.000  0.000
  Location.CTrt1  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 4Trt1  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 5Trt1  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 6Trt1 -0.656  0.000  0.000  0.000  0.000  0.000
Location ^ 7Trt1  0.000 -0.656  0.000  0.000  0.000  0.000
Location ^ 8Trt1  0.000  0.000 -0.656  0.000  0.000  0.000
  Location.LTrt2  0.000  0.000  0.000  0.000  0.000  0.000
  Location.QTrt2  0.000  0.000  0.000  0.000  0.000  0.000
  Location.CTrt2  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 4Trt2  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 5Trt2  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 6Trt2 -0.656  0.000  0.000  0.000  0.000  0.000
Location ^ 7Trt2  0.000 -0.656  0.000  0.000  0.000  0.000
Location ^ 8Trt2  0.000  0.000 -0.656  0.000  0.000  0.000
  Location.LTrt3  0.000  0.000  0.000  0.000  0.000  0.000
  Location.QTrt3  0.000  0.000  0.000  0.000  0.000  0.000
  Location.CTrt3  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 4Trt3  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 5Trt3  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 6Trt3 -0.656  0.000  0.000  0.000  0.000  0.000
Location ^ 7Trt3  0.000 -0.656  0.000  0.000  0.000  0.000
Location ^ 8Trt3  0.000  0.000 -0.656  0.000  0.000  0.000
                 LctLT1 LctQT1 LctCT1 Lc^4T1 Lc^5T1 Lc^6T1 
      Location.L                                          
      Location.Q                                          
      Location.C                                          
    Location ^ 4                                          
    Location ^ 5                                          
    Location ^ 6                                          
    Location ^ 7                                          
    Location ^ 8                                          
            Trt1                                          
            Trt2                                          
            Trt3                                          
  Location.LTrt1                                          
  Location.QTrt1  0.000                                   
  Location.CTrt1  0.000  0.000                            
Location ^ 4Trt1  0.000  0.000  0.000                     
Location ^ 5Trt1  0.000  0.000  0.000  0.000              
Location ^ 6Trt1  0.000  0.000  0.000  0.000  0.000       
Location ^ 7Trt1  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 8Trt1  0.000  0.000  0.000  0.000  0.000  0.000
  Location.LTrt2  0.500  0.000  0.000  0.000  0.000  0.000
  Location.QTrt2  0.000  0.500  0.000  0.000  0.000  0.000
  Location.CTrt2  0.000  0.000  0.500  0.000  0.000  0.000
Location ^ 4Trt2  0.000  0.000  0.000  0.500  0.000  0.000
Location ^ 5Trt2  0.000  0.000  0.000  0.000  0.500  0.000
Location ^ 6Trt2  0.000  0.000  0.000  0.000  0.000  0.500
Location ^ 7Trt2  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 8Trt2  0.000  0.000  0.000  0.000  0.000  0.000
  Location.LTrt3  0.500  0.000  0.000  0.000  0.000  0.000
  Location.QTrt3  0.000  0.500  0.000  0.000  0.000  0.000
  Location.CTrt3  0.000  0.000  0.500  0.000  0.000  0.000
Location ^ 4Trt3  0.000  0.000  0.000  0.500  0.000  0.000
Location ^ 5Trt3  0.000  0.000  0.000  0.000  0.500  0.000
Location ^ 6Trt3  0.000  0.000  0.000  0.000  0.000  0.500
Location ^ 7Trt3  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 8Trt3  0.000  0.000  0.000  0.000  0.000  0.000
                 Lc^7T1 Lc^8T1 LctLT2 LctQT2 LctCT2 Lc^4T2 
      Location.L                                          
      Location.Q                                          
      Location.C                                          
    Location ^ 4                                          
    Location ^ 5                                          
    Location ^ 6                                          
    Location ^ 7                                          
    Location ^ 8                                          
            Trt1                                          
            Trt2                                          
            Trt3                                          
  Location.LTrt1                                          
  Location.QTrt1                                          
  Location.CTrt1                                          
Location ^ 4Trt1                                          
Location ^ 5Trt1                                          
Location ^ 6Trt1                                          
Location ^ 7Trt1                                          
Location ^ 8Trt1  0.000                                   
  Location.LTrt2  0.000  0.000                            
  Location.QTrt2  0.000  0.000  0.000                     
  Location.CTrt2  0.000  0.000  0.000  0.000              
Location ^ 4Trt2  0.000  0.000  0.000  0.000  0.000       
Location ^ 5Trt2  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 6Trt2  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 7Trt2  0.500  0.000  0.000  0.000  0.000  0.000
Location ^ 8Trt2  0.000  0.500  0.000  0.000  0.000  0.000
  Location.LTrt3  0.000  0.000  0.500  0.000  0.000  0.000
  Location.QTrt3  0.000  0.000  0.000  0.500  0.000  0.000
  Location.CTrt3  0.000  0.000  0.000  0.000  0.500  0.000
Location ^ 4Trt3  0.000  0.000  0.000  0.000  0.000  0.500
Location ^ 5Trt3  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 6Trt3  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 7Trt3  0.500  0.000  0.000  0.000  0.000  0.000
Location ^ 8Trt3  0.000  0.500  0.000  0.000  0.000  0.000
                 Lc^5T2 Lc^6T2 Lc^7T2 Lc^8T2 LctLT3 LctQT3 
      Location.L                                          
      Location.Q                                          
      Location.C                                          
    Location ^ 4                                          
    Location ^ 5                                          
    Location ^ 6                                          
    Location ^ 7                                          
    Location ^ 8                                          
            Trt1                                          
            Trt2                                          
            Trt3                                          
  Location.LTrt1                                          
  Location.QTrt1                                          
  Location.CTrt1                                          
Location ^ 4Trt1                                          
Location ^ 5Trt1                                          
Location ^ 6Trt1                                          
Location ^ 7Trt1                                          
Location ^ 8Trt1                                          
  Location.LTrt2                                          
  Location.QTrt2                                          
  Location.CTrt2                                          
Location ^ 4Trt2                                          
Location ^ 5Trt2                                          
Location ^ 6Trt2  0.000                                   
Location ^ 7Trt2  0.000  0.000                            
Location ^ 8Trt2  0.000  0.000  0.000                     
  Location.LTrt3  0.000  0.000  0.000  0.000              
  Location.QTrt3  0.000  0.000  0.000  0.000  0.000       
  Location.CTrt3  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 4Trt3  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 5Trt3  0.500  0.000  0.000  0.000  0.000  0.000
Location ^ 6Trt3  0.000  0.500  0.000  0.000  0.000  0.000
Location ^ 7Trt3  0.000  0.000  0.500  0.000  0.000  0.000
Location ^ 8Trt3  0.000  0.000  0.000  0.500  0.000  0.000
                 LctCT3 Lc^4T3 Lc^5T3 Lc^6T3 Lc^7T3 
      Location.L                                   
      Location.Q                                   
      Location.C                                   
    Location ^ 4                                   
    Location ^ 5                                   
    Location ^ 6                                   
    Location ^ 7                                   
    Location ^ 8                                   
            Trt1                                   
            Trt2                                   
            Trt3                                   
  Location.LTrt1                                   
  Location.QTrt1                                   
  Location.CTrt1                                   
Location ^ 4Trt1                                   
Location ^ 5Trt1                                   
Location ^ 6Trt1                                   
Location ^ 7Trt1                                   
Location ^ 8Trt1                                   
  Location.LTrt2                                   
  Location.QTrt2                                   
  Location.CTrt2                                   
Location ^ 4Trt2                                   
Location ^ 5Trt2                                   
Location ^ 6Trt2                                   
Location ^ 7Trt2                                   
Location ^ 8Trt2                                   
  Location.LTrt3                                   
  Location.QTrt3                                   
  Location.CTrt3                                   
Location ^ 4Trt3  0.000                            
Location ^ 5Trt3  0.000  0.000                     
Location ^ 6Trt3  0.000  0.000  0.000              
Location ^ 7Trt3  0.000  0.000  0.000  0.000       
Location ^ 8Trt3  0.000  0.000  0.000  0.000  0.000

Standardized Within-Group Residuals:
     Min      Q1        Med      Q3    Max 
 -2.4871 -0.5169 -0.0009498 0.55424 3.2972

Number of Observations: 108
Number of Groups: 27 
S> fm2Mult <- update( fm1Mult, Adj ~ Location + Trt )
S> fm3Mult <- update( fm1Mult, Adj ~ Location )
S> fm4Mult <- update( fm1Mult, Adj ~ Trt )
S> fm5Mult <- update( fm1Mult, Adj ~ 1 )
S> anova( fm1Mult, fm2Mult, fm3Mult, fm5Mult )
        Model df     AIC    BIC  logLik   Test L.Ratio p-value 
fm1Mult     1 38 -11.146 90.775  43.573                       
fm2Mult     2 14 -24.504 13.046  26.252 1 vs 2  34.642  0.0739
fm3Mult     3 11  -0.956 28.547  11.478 2 vs 3  29.548  <.0001
fm5Mult     4  3  49.731 57.777 -21.865 3 vs 4  66.687  <.0001
S> anova( fm1Mult, fm2Mult, fm4Mult, fm5Mult )
        Model df     AIC    BIC  logLik   Test L.Ratio p-value 
fm1Mult     1 38 -11.146 90.775  43.573                       
fm2Mult     2 14 -24.504 13.046  26.252 1 vs 2  34.642  0.0739
fm4Mult     3  6  26.951 43.044  -7.476 2 vs 3  67.455  <.0001
fm5Mult     4  3  49.731 57.777 -21.865 3 vs 4  28.780  <.0001
S> ### AIC, BIC, and likelihood ratio tests all prefer model fm2Mult
S> summary( fm2Mult )
Linear mixed-effects model fit by maximum likelihood
 Data: Multilocation 
      AIC    BIC logLik 
  -24.504 13.046 26.252

Random effects:
 Formula:  ~ 1 | Grp
        (Intercept) Residual 
StdDev:    0.025894  0.18803

Fixed effects: Adj ~ Location + Trt 
               Value Std.Error DF t-value p-value 
 (Intercept)  2.8657  0.038744 78  73.964  <.0001
  Location.L  0.9176  0.059716 18  15.366  <.0001
  Location.Q  0.0026  0.059716 18   0.044  0.9655
  Location.C  0.0786  0.059716 18   1.316  0.2046
Location ^ 4 -0.1634  0.059716 18  -2.737  0.0135
Location ^ 5  0.0658  0.059716 18   1.102  0.2850
Location ^ 6  0.0471  0.059716 18   0.789  0.4406
Location ^ 7 -0.2422  0.059716 18  -4.055  0.0007
Location ^ 8  0.1235  0.059716 18   2.068  0.0534
        Trt1  0.0583  0.054280 78   1.075  0.2857
        Trt2 -0.1880  0.054280 78  -3.464  0.0009
        Trt3  0.0838  0.054280 78   1.544  0.1267
 Correlation: 
             (Intr) Lctn.L Lctn.Q Lctn.C Lctn^4 Lctn^5 Lctn^6 
  Location.L  0.0                                            
  Location.Q  0.0    0.0                                     
  Location.C  0.0    0.0    0.0                              
Location ^ 4  0.0    0.0    0.0    0.0                       
Location ^ 5  0.0    0.0    0.0    0.0    0.0                
Location ^ 6  0.0    0.0    0.0    0.0    0.0    0.0         
Location ^ 7  0.0    0.0    0.0    0.0    0.0    0.0    0.0  
Location ^ 8  0.0    0.0    0.0    0.0    0.0    0.0    0.0  
        Trt1 -0.7    0.0    0.0    0.0    0.0    0.0    0.0  
        Trt2 -0.7    0.0    0.0    0.0    0.0    0.0    0.0  
        Trt3 -0.7    0.0    0.0    0.0    0.0    0.0    0.0  
             Lctn^7 Lctn^8 Trt1 Trt2 
  Location.L                        
  Location.Q                        
  Location.C                        
Location ^ 4                        
Location ^ 5                        
Location ^ 6                        
Location ^ 7                        
Location ^ 8  0.0                   
        Trt1  0.0    0.0            
        Trt2  0.0    0.0    0.5     
        Trt3  0.0    0.0    0.5  0.5

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3    Max 
 -2.1253 -0.65599 -0.092189 0.61603 2.9345

Number of Observations: 108
Number of Groups: 27 
S> fm2RMult <- update( fm2Mult, method = "REML" ) # get REML estimates
S> summary( fm2RMult )
Linear mixed-effects model fit by REML
 Data: Multilocation 
     AIC    BIC logLik 
  24.196 60.097 1.9019

Random effects:
 Formula:  ~ 1 | Grp
        (Intercept) Residual 
StdDev:     0.07131  0.19161

Fixed effects: Adj ~ Location + Trt 
               Value Std.Error DF t-value p-value 
 (Intercept)  2.8657  0.039347 78  72.831  <.0001
  Location.L  0.9176  0.068954 18  13.307  <.0001
  Location.Q  0.0026  0.068954 18   0.038  0.9702
  Location.C  0.0786  0.068954 18   1.140  0.2693
Location ^ 4 -0.1634  0.068954 18  -2.370  0.0292
Location ^ 5  0.0658  0.068954 18   0.954  0.3526
Location ^ 6  0.0471  0.068954 18   0.683  0.5033
Location ^ 7 -0.2422  0.068954 18  -3.512  0.0025
Location ^ 8  0.1235  0.068954 18   1.791  0.0902
        Trt1  0.0583  0.052150 78   1.119  0.2667
        Trt2 -0.1880  0.052150 78  -3.605  0.0005
        Trt3  0.0838  0.052150 78   1.607  0.1122
 Correlation: 
             (Intr) Lctn.L Lctn.Q Lctn.C Lctn^4 Lctn^5 Lctn^6 
  Location.L  0.000                                          
  Location.Q  0.000  0.000                                   
  Location.C  0.000  0.000  0.000                            
Location ^ 4  0.000  0.000  0.000  0.000                     
Location ^ 5  0.000  0.000  0.000  0.000  0.000              
Location ^ 6  0.000  0.000  0.000  0.000  0.000  0.000       
Location ^ 7  0.000  0.000  0.000  0.000  0.000  0.000  0.000
Location ^ 8  0.000  0.000  0.000  0.000  0.000  0.000  0.000
        Trt1 -0.663  0.000  0.000  0.000  0.000  0.000  0.000
        Trt2 -0.663  0.000  0.000  0.000  0.000  0.000  0.000
        Trt3 -0.663  0.000  0.000  0.000  0.000  0.000  0.000
             Lctn^7 Lctn^8   Trt1   Trt2 
  Location.L                            
  Location.Q                            
  Location.C                            
Location ^ 4                            
Location ^ 5                            
Location ^ 6                            
Location ^ 7                            
Location ^ 8  0.000                     
        Trt1  0.000  0.000              
        Trt2  0.000  0.000  0.500       
        Trt3  0.000  0.000  0.500  0.500

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3    Max 
 -2.0603 -0.65887 -0.059717 0.60982 2.6211

Number of Observations: 108
Number of Groups: 27 
S> ### Treating the location as a random effect
S> fm1MultR <- lme( Adj ~ Trt, data = Multilocation, method = "ML",
+   random = list( Location = pdCompSymm( ~ Trt - 1 ), Block = ~ 1 ) )
S> summary( fm1MultR )
Linear mixed-effects model fit by maximum likelihood
 Data: Multilocation 
     AIC    BIC logLik 
  2.3564 23.813 6.8218

Random effects:
 Formula:  ~ Trt - 1 | Location
 Structure: Compound Symmetry
      StdDev  Corr             
Trt1 0.31933                  
Trt2 0.31933 0.992            
Trt3 0.31933 0.992 0.992      
Trt4 0.31933 0.992 0.992 0.992

 Formula:  ~ 1 | Block %in% Location
        (Intercept) Residual 
StdDev:    0.074962  0.18595

Fixed effects: Adj ~ Trt 
              Value Std.Error DF t-value p-value 
(Intercept)  2.8657   0.11538 78  24.837  <.0001
       Trt1  0.0583   0.05328 78   1.095  0.2769
       Trt2 -0.1880   0.05328 78  -3.529  0.0007
       Trt3  0.0838   0.05328 78   1.572  0.1199
 Correlation: 
     (Intr)   Trt1   Trt2 
Trt1 -0.231              
Trt2 -0.231  0.500       
Trt3 -0.231  0.500  0.500

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3   Max 
 -2.0977 -0.64759 -0.062211 0.58966 2.619

Number of Observations: 108
Number of Groups: 
 Location Block %in% Location 
        9                  27
S> fm2MultR <- update( fm1MultR, random = list( Location = ~ Trt - 1, Block = ~ 1 ))
Warning messages:
  ITERATION LIMIT REACHED WITHOUT OTHER CONVERGENCE in: ms( ~ \
	- logLik(lmeSt, lmePars), start = list(lmePars = c(coef( ....
S> anova( fm1MultR, fm2MultR )
         Model df    AIC    BIC logLik   Test L.Ratio p-value 
fm1MultR     1  8  2.356 23.813 6.8218                       
fm2MultR     2 16 15.825 58.739 8.0877 1 vs 2  2.5318  0.9602
S> ## No indication that a general variance-covariance is preferred to
S> ## a compound symmetry structure.
S> fm1RMultR <- update( fm1MultR, method = "REML" )
S> summary( fm1RMultR )
Linear mixed-effects model fit by REML
 Data: Multilocation 
     AIC    BIC   logLik 
  17.613 38.768 -0.80661

Random effects:
 Formula:  ~ Trt - 1 | Location
 Structure: Compound Symmetry
      StdDev Corr           
Trt1 0.34116               
Trt2 0.34116 0.98          
Trt3 0.34116 0.98 0.98     
Trt4 0.34116 0.98 0.98 0.98

 Formula:  ~ 1 | Block %in% Location
        (Intercept) Residual 
StdDev:    0.074962  0.18595

Fixed effects: Adj ~ Trt 
              Value Std.Error DF t-value p-value 
(Intercept)  2.8657   0.12009 78  23.863  <.0001
       Trt1  0.0583   0.05546 78   1.052  0.2960
       Trt2 -0.1880   0.05546 78  -3.390  0.0011
       Trt3  0.0838   0.05546 78   1.511  0.1349
 Correlation: 
     (Intr)   Trt1   Trt2 
Trt1 -0.231              
Trt2 -0.231  0.500       
Trt3 -0.231  0.500  0.500

Standardized Within-Group Residuals:
     Min       Q1       Med     Q3    Max 
 -1.9837 -0.67435 -0.030874 0.5388 2.6272

Number of Observations: 108
Number of Groups: 
 Location Block %in% Location 
        9                  27
S> c( 0.34116, 0.07497, 0.18596)^2  # compare with estimates, p. 84
[1] 0.1163901 0.0056205 0.0345811
