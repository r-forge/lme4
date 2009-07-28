library(mlmRev)
options(show.signif.stars = FALSE)
lmer(mAch ~ meanses*cses + sector*cses + (cses|school), Hsb82)
q("no")
