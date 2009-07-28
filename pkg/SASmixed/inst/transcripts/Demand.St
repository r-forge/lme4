S> ### $Id: Demand.St,v 1.1 1999/10/13 00:50:09 saikat Exp $
S> ### Analysis of the per capita demand deposits data given as data
S> ### set 3.6 in "SAS System for Mixed Models"
S> names( Demand )
[1] "State" "Year"  "d"     "y"     "rd"    "rt"    "rs"   
[8] "grp"  
S> formula( Demand )
d ~ Year | State
S> # the "grp" factor is a dummy factor with only one level.
S> unique( Demand$grp )
[1] 1
S> # Crossed random-effects factors have to be created by pdIdent applied to
S> # the indicator variables and joined together by pdBlocked.
S> fm1Demand <- lme( log(d) ~ log(y) + log(rd) + log(rt) + log(rs), data = Demand,
+   random = list(grp = pdBlocked(list(pdIdent(~ State - 1), pdIdent(~ Year - 1)))),
+   method = "REML" )
S> summary( fm1Demand )        # compare to output 3.13, p. 132
Linear mixed-effects model fit by REML
 Data: Demand 
      AIC     BIC logLik 
  -224.17 -205.95 120.08

Random effects:
 Composite Structure: Blocked

 Block 1: StateWA, StateFL, StateCA, StateTX, StateIL, StateDC, StateNY
 Formula:  ~ State - 1 | grp
 Structure: Multiple of an Identity
        StateWA StateFL StateCA StateTX StateIL StateDC StateNY 
StdDev: 0.17177 0.17177 0.17177 0.17177 0.17177 0.17177 0.17177

 Block 2: Year1949, Year1950, Year1951, Year1952, Year1953, Year1954, Year1955, Year1956, Year1957, Year1958, Year1959
 Formula:  ~ Year - 1 | grp
 Structure: Multiple of an Identity
        Year1949 Year1950 Year1951 Year1952 Year1953 Year1954 
StdDev: 0.016267 0.016267 0.016267 0.016267 0.016267 0.016267
        Year1955 Year1956 Year1957 Year1958 Year1959 Residual 
StdDev: 0.016267 0.016267 0.016267 0.016267 0.016267 0.033421

Fixed effects: log(d) ~ log(y) + log(rd) + log(rt) + log(rs) 
              Value Std.Error DF t-value p-value 
(Intercept) -1.2839   0.72343 72  -1.775  0.0802
     log(y)  1.0698   0.10393 72  10.294  <.0001
    log(rd) -0.2953   0.05246 72  -5.629  <.0001
    log(rt)  0.0399   0.02789 72   1.430  0.1570
    log(rs) -0.3267   0.11438 72  -2.857  0.0056
 Correlation: 
        (Intr) log(y) lg(rd) lg(rt) 
 log(y) -0.976                     
log(rd)  0.383 -0.227              
log(rt)  0.077 -0.062 -0.337       
log(rs)  0.444 -0.600 -0.270 -0.323

Standardized Within-Group Residuals:
     Min       Q1      Med      Q3    Max 
 -2.1382 -0.59488 0.045884 0.47185 2.5816

Number of Observations: 77
Number of Groups: 1 
