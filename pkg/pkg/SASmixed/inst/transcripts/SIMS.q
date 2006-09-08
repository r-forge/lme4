### $Id: SIMS.q,v 1.1 1999/10/13 00:50:09 saikat Exp $
### Analysis of the Second International Mathematics Study (SIMS)
### described in section 7.2 of "SAS System for Mixed Models"
plot(SIMS, layout = c(10, 10, 2) )  # high level of variability within classes
## to record the execution time for comparison the fit was done as
## unix.time(assign("fm1RSIMS", lme(Gain ~ Pretot, SIMS, ~ Pretot | Class, REML = TRUE, control = list(msVerbose = TRUE))))
## Since unix.time is not available on Windows, we give the equivalent
## statement here
fm1RSIMS <- lme(Gain ~ Pretot, data = SIMS, 
    random = ~ Pretot | Class, control = list(msVerbose = TRUE))
## Timing done on an UltraSPARC 2 running S-PLUS 3.4 for Solaris.
## Your mileage may vary.
summary(fm1RSIMS)              # compare to output 7.4, p. 262
intervals( fm1RSIMS )
