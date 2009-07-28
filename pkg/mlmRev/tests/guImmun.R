library(mlmRev)
options(show.signif.stars = FALSE)
(fm <- lmer(immun ~ kid2p + mom25p + ord + ethn + momEd +
            husEd + momWork + rural + pcInd81 + (1|mom) + (1|comm),
            guImmun, family = binomial))
q("no")
