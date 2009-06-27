###################################################
### chunk number 1: preliminaries
###################################################
options(width=65,show.signif.stars=FALSE)
library(lattice)
library(Matrix)
library(lme4)
data(Machines, package = "MEMSS")
Machines <- Machines[with(Machines,order(Worker, Machine)),]
lattice.options(default.theme = function() standard.theme())
#lattice.options(default.theme = function() standard.theme(color=FALSE))
if (exists("ratbrain.rda")) {
    load("ratbrain.rda")
} else {
    ratbrain <- within(read.delim("http://www-personal.umich.edu/~bwest/rat_brain.dat"),
                   {
                       treatment <- factor(treatment, 
                                           labels = c("Basal", "Carbachol"))
                       region <- factor(region, levels = c(3, 1, 2),
                                        labels = c("VDB", "BST", "LS"))
                   })
    save(ratbrain, file = "ratbrain.rda")
}


###################################################
### chunk number 2: Machinesplot
###################################################
print(dotplot(reorder(Worker, score) ~ score, Machines,
              groups = Machine, 
              xlab = "Quality and productivity score", 
              ylab = "Worker", type = c("p","a"), 
              auto.key = list(columns = 3, lines = TRUE), 
              jitter.y = TRUE))


###################################################
### chunk number 3: MachinesModel1
###################################################
print(fm1 <- lmer(score ~ Machine + (1|Worker) + (1|Worker:Machine),
                  Machines), corr = FALSE)


###################################################
### chunk number 4: fm1Z
###################################################
print(image(fm1@Zt, sub = NULL))


###################################################
### chunk number 5: MachinesModel1
###################################################
print(fm2 <- lmer(score ~ Machine + (0+Machine|Worker), Machines),
      corr = FALSE)


###################################################
### chunk number 6: fm2Z
###################################################
print(image(fm2@Zt, sub = NULL, xlab = NULL, ylab = NULL))


###################################################
### chunk number 7: Machinesanova
###################################################
fm2M <- update(fm2, REML = FALSE)
fm1M <- update(fm1, REML = FALSE)
anova(fm2M, fm1M)


###################################################
### chunk number 8: MachineswoB6
###################################################
Machines1 <- subset(Machines, !(Worker == "6" & Machine == 'B'))
xtabs(~ Machine + Worker, Machines1)


###################################################
### chunk number 9: Machines1plot
###################################################
print(dotplot(reorder(Worker, score) ~ score, Machines1,
              groups = Machine, 
              xlab = "Quality and productivity score", 
              ylab = "Worker", type = c("p","a"), 
              auto.key = list(columns = 3, lines = TRUE), 
              jitter.y = TRUE))


###################################################
### chunk number 10: fm1aM
###################################################
fm1aM <- lmer(score ~ Machine + (1|Worker) + (1|Worker:Machine), Machines1, REML = FALSE)
fm2aM <- lmer(score ~ Machine + (0 + Machine|Worker), Machines1, REML = FALSE)
anova(fm2aM, fm1aM)


###################################################
### chunk number 11: ratbraindot1
###################################################
print(dotplot(factor(region, levels = c("BST", "LS", "VDB")) ~ activate|treatment,
              ratbrain, groups = animal, type = c("p","a"),
              strip = FALSE, strip.left = TRUE,
              xlab = "Activation (mean optical density)",
              ylab = "Region of brain",
              layout = c(1,2),
              auto.key = list(columns = 5, lines = TRUE, points = FALSE)))
              


###################################################
### chunk number 12: ratbraindot
###################################################
print(dotplot(region ~ activate|animal,
              ratbrain, groups = treatment, type = c("p","a"),
              strip = FALSE, strip.left = TRUE,
              xlab = "Activation (mean optical density)",
              ylab = "Region",
              layout = c(1,5),
              auto.key = list(columns = 2, lines = TRUE, points = FALSE)))


###################################################
### chunk number 13: ratbraindat
###################################################
str(ratbrain)


###################################################
### chunk number 14: m51
###################################################
print(m51 <- lmer(activate ~ region * treatment + (1|animal), ratbrain), corr = FALSE)


###################################################
### chunk number 15: m52
###################################################
print(m52 <- lmer(activate ~ region * treatment + (treatment|animal), ratbrain), corr = FALSE)


###################################################
### chunk number 16: m52
###################################################
print(m52a <- lmer(activate ~ region * treatment + (0+treatment|animal), ratbrain), corr = FALSE)


###################################################
### chunk number 17: m54
###################################################
print(m54 <- lmer(activate ~ region * treatment + (1|animal) +
                  (1|animal:treatment), ratbrain), corr = FALSE)


###################################################
### chunk number 18: predintact
###################################################
print(dotplot(ranef(m51, post = TRUE),
              strip = FALSE, strip.left = TRUE)$animal,
      split = c(1,1,1,3), more = TRUE)
print(dotplot(ranef(m52, post = TRUE),
              strip = FALSE, strip.left = TRUE,
              scales = list(x = list(relation = "free")))$animal,
      split = c(1,2,1,3), more = TRUE)
print(dotplot(ranef(m52a, post = TRUE),
              strip = FALSE, strip.left = TRUE,
              scales = list(x = list(relation = "free")))$animal,
      split = c(1,3,1,3))


###################################################
### chunk number 19: anovabrain
###################################################
m51M <- update(m51, REML = FALSE)
m52M <- update(m52, REML = FALSE)
anova(m51M, m52M)


###################################################
### chunk number 20: ftable
###################################################
ftable(xtabs(~ treatment + region + animal, ratbrain))


###################################################
### chunk number 21: trt
###################################################
ftable(atab <- xtabs(activate ~ treatment + animal + region, ratbrain))


###################################################
### chunk number 22: stratab
###################################################
str(atab)


###################################################
### chunk number 23: diffstab
###################################################
(diffs <- as.table(apply(atab, 2:3, diff)))


###################################################
### chunk number 24: diffsfr
###################################################
str(diffs <- as.data.frame(diffs))
names(diffs)[3] <- "actdiff"


###################################################
### chunk number 25: actdiffdat
###################################################
print(dotplot(reorder(animal, actdiff) ~ actdiff, diffs, groups = region,
              xlab = "Difference in activation with Carbachol from Basal state",
              ylab = "Animal", type = c("p","a"),
              auto.key = list(columns = 3, lines = TRUE, points = FALSE)))


###################################################
### chunk number 26: dm1
###################################################
print(dm1 <- lmer(actdiff ~ region + (1|animal), diffs))


###################################################
### chunk number 27: aovmod
###################################################
summary(m52f <- aov(activate ~ animal * treatment + region * treatment,
                    ratbrain))
anova(m52)


