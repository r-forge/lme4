###################################################
### chunk number 1: initial
###################################################
options(width=80, show.signif.stars = FALSE)
library(lattice)
lattice.options(default.theme = standard.theme(color = FALSE))
set.seed(123454321)


###################################################
### chunk number 2: loadclassroom
###################################################
load("classroom.rda")
library(lattice)


###################################################
### chunk number 3: histogram
###################################################
print(histogram(~ mathkind, classroom,
                xlab = "Kindergarten mathematics score"))


###################################################
### chunk number 4: density
###################################################
print(densityplot(~ mathkind, classroom,
                  xlab = "Kindergarten mathematics score"))


###################################################
### chunk number 5: densitysex
###################################################
print(densityplot(~ mathkind, classroom, groups = sex,
                  auto.key = list(columns = 2), plot.points = FALSE,
                  xlab = "Kindergarten mathematics score"))


###################################################
### chunk number 6: densityminority
###################################################
print(densityplot(~ mathkind, classroom, groups = minority,
                  auto.key = list(columns = 2), plot.points = FALSE,
                  xlab = "Kindergarten mathematics score"))


###################################################
### chunk number 7: xyplot
###################################################
print(xyplot(mathgain ~ mathkind, classroom,
             type = c("g","p","smooth"),
             ylab = "Change in mathematics score",
             xlab = "Kindergarten mathematics score"),
      split = c(1,1,2,1), more = TRUE)
print(xyplot(mathgain ~ mathkind, classroom,
             type = c("g","p","r"),
             ylab = "Change in mathematics score",
             xlab = "Kindergarten mathematics score"),
      split = c(2,1,2,1))


###################################################
### chunk number 8: corr
###################################################
with(classroom, cor(mathkind, mathgain))


###################################################
### chunk number 9: math1
###################################################
classroom <- within(classroom, math1 <- mathkind + mathgain)
with(classroom, cor(mathkind, math1))


###################################################
### chunk number 10: xyplot1
###################################################
print(xyplot(math1 ~ mathkind, classroom,
             type = c("g","p","smooth"), aspect = "iso",
             ylab = "Grade 1 mathematics score",
             xlab = "Kindergarten mathematics score"),
      split = c(1,1,2,1), more = TRUE)
print(xyplot(math1 ~ mathkind, classroom,
             type = c("g","p","r"), aspect = "iso",
             ylab = "Grade 1 mathematics score",
             xlab = "Kindergarten mathematics score"),
      split = c(2,1,2,1))


###################################################
### chunk number 11: xyplot2
###################################################
print(xyplot(math1 ~ mathkind|sex, classroom,
             type = c("g","p","r"), aspect = "iso",
             ylab = "Grade 1 mathematics score",
             xlab = "Kindergarten mathematics score"))


###################################################
### chunk number 12: xyplot3
###################################################
print(xyplot(math1 ~ mathkind|minority, classroom,
             type = c("g","p","r"), aspect = "iso",
             ylab = "Grade 1 mathematics score",
             xlab = "Kindergarten mathematics score"))


###################################################
### chunk number 13: xyplot4
###################################################
print(xyplot(math1 ~ mathkind|minority*sex, classroom,
             type = c("g","p","r"), aspect = "iso",
             ylab = "Grade 1 mathematics score",
             xlab = "Kindergarten mathematics score"))


###################################################
### chunk number 14: xtabs
###################################################
xtabs(~ classid, classroom, schoolid == 11, drop = TRUE)


###################################################
### chunk number 15: dotplot1
###################################################
print(dotplot(classid ~ mathgain, classroom, subset = schoolid == 11,
              pch = 21, jitter.y = TRUE,
              xlab = "Gain in mathematics score",
              ylab = "Grade 1 classroom (school 11 only)"))


###################################################
### chunk number 16: dotplot2
###################################################
print(dotplot(reorder(classid,mathgain) ~ mathgain, classroom,
              subset = schoolid == 11, type = c("p","a"),
              pch = 21, jitter.y = TRUE,
              xlab = "Gain in mathematics score",
              ylab = "Grade 1 classroom (school 11 only)"))


