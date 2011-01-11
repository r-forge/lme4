###################################################
### chunk number 1: initial
###################################################
#line 20 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
options(width=80, show.signif.stars = FALSE)
library(lattice)
lattice.options(default.theme = standard.theme(color = FALSE))
set.seed(123454321)


###################################################
### chunk number 2: loadclassroom
###################################################
#line 29 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
load("classroom.rda")
library(lattice)


###################################################
### chunk number 3: histogram
###################################################
#line 37 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
print(histogram(~ mathkind, classroom,
                xlab = "Kindergarten mathematics score"))


###################################################
### chunk number 4: density
###################################################
#line 44 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
print(densityplot(~ mathkind, classroom,
                  xlab = "Kindergarten mathematics score"))


###################################################
### chunk number 5: densitysex
###################################################
#line 55 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
print(densityplot(~ mathkind, classroom, groups = sex,
                  auto.key = list(columns = 2), plot.points = FALSE,
                  xlab = "Kindergarten mathematics score"))


###################################################
### chunk number 6: densityminority
###################################################
#line 65 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
print(densityplot(~ mathkind, classroom, groups = minority,
                  auto.key = list(columns = 2), plot.points = FALSE,
                  xlab = "Kindergarten mathematics score"))


###################################################
### chunk number 7: xyplot
###################################################
#line 79 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
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
#line 94 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
with(classroom, cor(mathkind, mathgain))


###################################################
### chunk number 9: math1
###################################################
#line 100 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
classroom <- within(classroom, math1 <- mathkind + mathgain)
with(classroom, cor(mathkind, math1))


###################################################
### chunk number 10: xyplot1
###################################################
#line 108 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
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
#line 125 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
print(xyplot(math1 ~ mathkind|sex, classroom,
             type = c("g","p","r"), aspect = "iso",
             ylab = "Grade 1 mathematics score",
             xlab = "Kindergarten mathematics score"))


###################################################
### chunk number 12: xyplot3
###################################################
#line 136 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
print(xyplot(math1 ~ mathkind|minority, classroom,
             type = c("g","p","r"), aspect = "iso",
             ylab = "Grade 1 mathematics score",
             xlab = "Kindergarten mathematics score"))


###################################################
### chunk number 13: xyplot4
###################################################
#line 146 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
print(xyplot(math1 ~ mathkind|minority*sex, classroom,
             type = c("g","p","r"), aspect = "iso",
             ylab = "Grade 1 mathematics score",
             xlab = "Kindergarten mathematics score"))


###################################################
### chunk number 14: xtabs
###################################################
#line 155 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
xtabs(~ classid, classroom, schoolid == 11, drop = TRUE)


###################################################
### chunk number 15: dotplot1
###################################################
#line 161 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
print(dotplot(classid ~ mathgain, classroom, subset = schoolid == 11,
              pch = 21, jitter.y = TRUE,
              xlab = "Gain in mathematics score",
              ylab = "Grade 1 classroom (school 11 only)"))


###################################################
### chunk number 16: dotplot2
###################################################
#line 171 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises2.Rnw"
print(dotplot(reorder(classid,mathgain) ~ mathgain, classroom,
              subset = schoolid == 11, type = c("p","a"),
              pch = 21, jitter.y = TRUE,
              xlab = "Gain in mathematics score",
              ylab = "Grade 1 classroom (school 11 only)"))


