S> ### $Id: PBIB.St,v 1.1 1999/10/13 00:50:09 saikat Exp $
S> ### Analysis of the partially balanced incomplete blocked PBIB data
S> options( contrasts = c(factor = "contr.SAS", ordered = "contr.poly") )
S> formula( PBIB )
response ~ Treatment | Block
S> names( PBIB )
[1] "response"  "Treatment" "Block"    
S> sapply( PBIB, data.class )
  response Treatment     Block 
 "numeric" "factor"  "ordered"
S> fm1PBIB <- lme( response ~ Treatment, data = PBIB, random = ~ 1 | Block,
+                method = "ML")
S> summary( fm1PBIB )
Linear mixed-effects model fit by maximum likelihood
 Data: PBIB 
     AIC    BIC  logLik 
  56.571 92.174 -11.285

Random effects:
 Formula:  ~ 1 | Block
        (Intercept) Residual 
StdDev:      0.2122  0.24571

Fixed effects: response ~ Treatment 
              Value Std.Error DF t-value p-value 
(Intercept)  2.8573   0.16577 31  17.236  <.0001
 Treatment1 -0.0346   0.21700 31  -0.159  0.8744
 Treatment2 -0.4507   0.21700 31  -2.077  0.0462
 Treatment3 -0.4026   0.21700 31  -1.855  0.0731
 Treatment4 -0.0923   0.21700 31  -0.425  0.6735
 Treatment5 -0.0520   0.22253 31  -0.234  0.8169
 Treatment6  0.0521   0.21700 31   0.240  0.8119
 Treatment7 -0.0718   0.21700 31  -0.331  0.7429
 Treatment8 -0.0700   0.21700 31  -0.323  0.7491
 Treatment9  0.0419   0.21700 31   0.193  0.8483
Treatment10 -0.3616   0.22253 31  -1.625  0.1143
Treatment11  0.0330   0.21700 31   0.152  0.8803
Treatment12  0.1879   0.21700 31   0.866  0.3933
Treatment13 -0.2297   0.21700 31  -1.059  0.2980
Treatment14 -0.3581   0.21700 31  -1.650  0.1090
 Correlation: 
            (Intr) Trtmn1 Trtmn2 Trtmn3 Trtmn4 Trtmn5 Trtmn6 
 Treatment1 -0.655                                          
 Treatment2 -0.655  0.500                                   
 Treatment3 -0.655  0.500  0.500                            
 Treatment4 -0.655  0.500  0.500  0.500                     
 Treatment5 -0.671  0.513  0.513  0.513  0.513              
 Treatment6 -0.655  0.474  0.500  0.500  0.500  0.513       
 Treatment7 -0.655  0.500  0.474  0.500  0.500  0.513  0.500
 Treatment8 -0.655  0.500  0.500  0.474  0.500  0.513  0.500
 Treatment9 -0.655  0.500  0.500  0.500  0.474  0.513  0.500
Treatment10 -0.671  0.513  0.513  0.513  0.513  0.500  0.513
Treatment11 -0.655  0.474  0.500  0.500  0.500  0.513  0.474
Treatment12 -0.655  0.500  0.474  0.500  0.500  0.513  0.500
Treatment13 -0.655  0.500  0.500  0.474  0.500  0.513  0.500
Treatment14 -0.655  0.500  0.500  0.500  0.474  0.513  0.500
            Trtmn7 Trtmn8 Trtmn9 Trtm10 Trtm11 Trtm12 Trtm13 
 Treatment1                                                 
 Treatment2                                                 
 Treatment3                                                 
 Treatment4                                                 
 Treatment5                                                 
 Treatment6                                                 
 Treatment7                                                 
 Treatment8  0.500                                          
 Treatment9  0.500  0.500                                   
Treatment10  0.513  0.513  0.513                            
Treatment11  0.500  0.500  0.500  0.513                     
Treatment12  0.474  0.500  0.500  0.513  0.500              
Treatment13  0.500  0.474  0.500  0.513  0.500  0.500       
Treatment14  0.500  0.500  0.474  0.513  0.500  0.500  0.500

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3    Max 
 -1.7419 -0.56754 0.0096936 0.44333 3.1061

Number of Observations: 60
Number of Groups: 15 
S> plot( fm1PBIB, resid(.) ~ fitted(.) | Block )
S> plot( fm1PBIB, resid(.) ~ fitted(.) | Treatment, inner = ~ Block )
Error in panel.xyplot(x, y, inner = struct..: Can't pass data of
	mode "call" to old-S function
Dumped
S> fm1RPBIB <- update( fm1PBIB, method = "REML" )
S> summary( fm1RPBIB )    # compare with output 1.7  pp. 24-25
Linear mixed-effects model fit by REML
 Data: PBIB 
     AIC   BIC  logLik 
  85.985 116.7 -25.992

Random effects:
 Formula:  ~ 1 | Block
        (Intercept) Residual 
StdDev:     0.21569   0.2925

Fixed effects: response ~ Treatment 
              Value Std.Error DF t-value p-value 
(Intercept)  2.8592   0.16641 31  17.182  <.0001
 Treatment1 -0.0417   0.22206 31  -0.188  0.8522
 Treatment2 -0.4539   0.22206 31  -2.044  0.0495
 Treatment3 -0.4043   0.22206 31  -1.821  0.0783
 Treatment4 -0.0754   0.22206 31  -0.340  0.7365
 Treatment5 -0.0543   0.22720 31  -0.239  0.8126
 Treatment6  0.0515   0.22206 31   0.232  0.8183
 Treatment7 -0.0702   0.22206 31  -0.316  0.7539
 Treatment8 -0.0776   0.22206 31  -0.350  0.7290
 Treatment9  0.0321   0.22206 31   0.144  0.8861
Treatment10 -0.3682   0.22720 31  -1.620  0.1153
Treatment11  0.0395   0.22206 31   0.178  0.8601
Treatment12  0.1936   0.22206 31   0.872  0.3900
Treatment13 -0.2415   0.22206 31  -1.087  0.2853
Treatment14 -0.3679   0.22206 31  -1.657  0.1076
 Correlation: 
            (Intr) Trtmn1 Trtmn2 Trtmn3 Trtmn4 Trtmn5 Trtmn6 
 Treatment1 -0.667                                          
 Treatment2 -0.667  0.500                                   
 Treatment3 -0.667  0.500  0.500                            
 Treatment4 -0.667  0.500  0.500  0.500                     
 Treatment5 -0.683  0.512  0.512  0.512  0.512              
 Treatment6 -0.667  0.477  0.500  0.500  0.500  0.512       
 Treatment7 -0.667  0.500  0.477  0.500  0.500  0.512  0.500
 Treatment8 -0.667  0.500  0.500  0.477  0.500  0.512  0.500
 Treatment9 -0.667  0.500  0.500  0.500  0.477  0.512  0.500
Treatment10 -0.683  0.512  0.512  0.512  0.512  0.500  0.512
Treatment11 -0.667  0.477  0.500  0.500  0.500  0.512  0.477
Treatment12 -0.667  0.500  0.477  0.500  0.500  0.512  0.500
Treatment13 -0.667  0.500  0.500  0.477  0.500  0.512  0.500
Treatment14 -0.667  0.500  0.500  0.500  0.477  0.512  0.500
            Trtmn7 Trtmn8 Trtmn9 Trtm10 Trtm11 Trtm12 Trtm13 
 Treatment1                                                 
 Treatment2                                                 
 Treatment3                                                 
 Treatment4                                                 
 Treatment5                                                 
 Treatment6                                                 
 Treatment7                                                 
 Treatment8  0.500                                          
 Treatment9  0.500  0.500                                   
Treatment10  0.512  0.512  0.512                            
Treatment11  0.500  0.500  0.500  0.512                     
Treatment12  0.477  0.500  0.500  0.512  0.500              
Treatment13  0.500  0.477  0.500  0.512  0.500  0.500       
Treatment14  0.500  0.500  0.477  0.512  0.500  0.500  0.500

Standardized Within-Group Residuals:
     Min       Q1       Med      Q3    Max 
 -1.5247 -0.45435 -0.013268 0.42568 2.7443

Number of Observations: 60
Number of Groups: 15 
S> ## Testing for significant fixed effects for Treatment involved re-fitting
S> ## the ML fit and comparing with anova().
S> fm2PBIB <- update( fm1PBIB, response ~ 1 )
S> anova( fm1PBIB, fm2PBIB )     
        Model df    AIC    BIC  logLik   Test L.Ratio p-value 
fm1PBIB     1 17 56.571 92.174 -11.285                       
fm2PBIB     2  3 52.152 58.435 -23.076 1 vs 2  23.581  0.0514
S> ## The p-value is considerably smaller than that for the F test in PROC MIXED
