S> ### $Id: Genetics.St,v 1.1 1999/10/13 00:50:09 saikat Exp $
S> ### Analysis of the heritability data given as data set 4.5 in
S> ### "SAS System for Mixed Models"
S> options(contrasts = c(factor = "contr.SAS", ordered = "contr.poly"))
S> formula( Genetics )
Yield ~ 1 | Location/Block
S> plot( Genetics )
S> fm1Gen <- lme( Yield ~ 1, data = Genetics, method = "ML", 
+    random = list(Location = pdCompSymm(~ Family - 1), Block = ~ 1) )
S> summary( fm1Gen )
Linear mixed-effects model fit by maximum likelihood
 Data: Genetics 
     AIC    BIC  logLik 
  499.71 510.19 -244.86

Random effects:
 Formula:  ~ Family - 1 | Location
 Structure: Compound Symmetry
        StdDev Corr                
Family1 25.951                    
Family2 25.951 0.61               
Family3 25.951 0.61 0.61          
Family4 25.951 0.61 0.61 0.61     
Family5 25.951 0.61 0.61 0.61 0.61

 Formula:  ~ 1 | Block %in% Location
        (Intercept) Residual 
StdDev:       9.451   7.2004

Fixed effects: Yield ~ 1 
             Value Std.Error DF t-value p-value 
(Intercept) 209.13    11.234 48  18.616  <.0001

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3    Max 
 -1.7827 -0.50294 -0.073372 0.62515 1.3794

Number of Observations: 60
Number of Groups: 
 Location Block %in% Location 
        4                  12
S> summary( update( fm1Gen, method = "REML" ) )
Linear mixed-effects model fit by REML
 Data: Genetics 
     AIC    BIC  logLik 
  492.92 503.31 -241.46

Random effects:
 Formula:  ~ Family - 1 | Location
 Structure: Compound Symmetry
        StdDev  Corr                   
Family1 28.964                        
Family2 28.964 0.687                  
Family3 28.964 0.687 0.687            
Family4 28.964 0.687 0.687 0.687      
Family5 28.964 0.687 0.687 0.687 0.687

 Formula:  ~ 1 | Block %in% Location
        (Intercept) Residual 
StdDev:       9.451   7.2004

Fixed effects: Yield ~ 1 
             Value Std.Error DF t-value p-value 
(Intercept) 209.13    12.863 48  16.258  <.0001

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3    Max 
 -1.7806 -0.51069 -0.070189 0.61932 1.3822

Number of Observations: 60
Number of Groups: 
 Location Block %in% Location 
        4                  12
