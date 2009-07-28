### $Id: Genetics.q,v 1.1 1999/10/13 00:50:09 saikat Exp $
### Analysis of the heritability data given as data set 4.5 in
### "SAS System for Mixed Models"
options(contrasts = c(factor = "contr.SAS", ordered = "contr.poly"))
formula( Genetics )
plot( Genetics )
fm1Gen <- lme( Yield ~ 1, data = Genetics, method = "ML", 
   random = list(Location = pdCompSymm(~ Family - 1), Block = ~ 1) )
summary( fm1Gen )
summary( update( fm1Gen, method = "REML" ) )
