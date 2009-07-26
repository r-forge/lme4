library(lattice)
load("classroom.rda")
print(histogram(~ mathkind, classroom,
                xlab = "Kindergarten mathematics score"))
print(densityplot(~ mathkind, classroom,
                  xlab = "Kindergarten mathematics score"))
print(densityplot(~ mathkind, classroom, groups = sex,
                  auto.key = list(columns = 2), plot.points = FALSE,
                  xlab = "Kindergarten mathematics score"))
print(densityplot(~ mathkind, classroom, groups = minority,
                  auto.key = list(columns = 2), plot.points = FALSE,
                  xlab = "Kindergarten mathematics score"))
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
with(classroom, cor(mathkind, mathgain))
classroom <- within(classroom, math1 <- mathkind + mathgain)
with(classroom, cor(mathkind, math1))
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
print(xyplot(math1 ~ mathkind|sex, classroom,
             type = c("g","p","r"), aspect = "iso",
             ylab = "Grade 1 mathematics score",
             xlab = "Kindergarten mathematics score"))
print(xyplot(math1 ~ mathkind|minority, classroom,
             type = c("g","p","r"), aspect = "iso",
             ylab = "Grade 1 mathematics score",
             xlab = "Kindergarten mathematics score"))
print(xyplot(math1 ~ mathkind|minority*sex, classroom,
             type = c("g","p","r"), aspect = "iso",
             ylab = "Grade 1 mathematics score",
             xlab = "Kindergarten mathematics score"))
xtabs(~ classid, classroom, schoolid == 11, drop = TRUE)
