library(mlmRev)
options(show.signif.stars = FALSE)
Early$tos <- Early$age - 0.5
(fm1 <- lmer(cog ~ tos * trt + (tos|id), Early))
q("no")
