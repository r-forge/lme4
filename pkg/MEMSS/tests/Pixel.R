library(MEMSS)
options(show.signif.stars = FALSE)
m1 <- lmer(pixel ~ day + I(day^2) + (1|Dog:Side) + (day|Dog),
           Pixel, verbose = TRUE)
print(m1, corr = FALSE)
ranef(m1)
q("no")
