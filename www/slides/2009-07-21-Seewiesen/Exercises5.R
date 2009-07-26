###################################################
### chunk number 1: initial
###################################################
options(width=80, show.signif.stars = FALSE)
library(MEMSS)
lattice.options(default.theme = standard.theme(color = FALSE))
set.seed(123454321)


###################################################
### chunk number 2: orthofem
###################################################
print(xyplot(distance ~ age|Subject, Orthodont, subset = Sex == "Female",
             index.cond = function(x,y) y[x == 8],
             aspect = 'xy', layout = c(11,1), type = c("g","p","r"),
             xlab = "Age (yr)",
             ylab = "Distance from pituitary to pterygomaxillary fissure (mm)"))

str(Orthodont)
Orth <- within(Orthodont, agec <- age - 11)
(fm1 <- lmer(distance ~ 1 + agec + (1 + agec|Subject), Orth, subset = Sex == "Female"))
dotplot(ranef(fm1, post = TRUE), scales = list(x = list(relation = "free")))
(fm2 <- lmer(distance ~ agec * Sex + (1 + agec|Subject), Orth))
