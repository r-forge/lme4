library(mlmRev)
options(show.signif.stars = FALSE)
cntr <- list()
#cntr <- list(EMverbose = TRUE, msVerbose = 1) # uncomment for verbose output
(fm1 <- lmer(attain ~ verbal * sex + (1|primary) + (1|second), ScotsSec,
             control = cntr))
(fm2 <- lmer(attain ~ verbal * sex + (1|primary) + (1|second), ScotsSec,
             control = c(cntr, list(niterEM = 40))))
## fm1 and fm2 should be essentially identical when optimizing with nlminb
## The fits are substantially different when optimizing with optim
q("no")
