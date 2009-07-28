library(mlmRev)
options(show.signif.stars = FALSE)
(fm1 <- lmer(langPOST ~ IQ.ver.cen+avg.IQ.ver.cen+(IQ.ver.cen|schoolNR),
             bdf))
q("no")

