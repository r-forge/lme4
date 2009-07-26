library(MEMSS)
str(Rail)
dotplot(reorder(Rail, travel) ~ travel, Rail, xlab = "Travel time",
        type = c("p", "a"), pch = 21, aspect = 0.3)
(fm1 <- lmer(travel ~ 1 + (1|Rail), Rail))
dotplot(ranef(fm1, post = TRUE), aspect = 0.25)
(fm1M <- update(fm1, REML = FALSE))
dotplot(ranef(fm1M, post = TRUE), aspect = 0.25)
rbind(ranef(fm1, drop = TRUE)[[1]], ranef(fm1M, drop = TRUE)[[1]])

dotplot(reorder(Subject, effort) ~ effort, ergoStool, groups = Type,
        type = c("p", "a"), auto.key = list(columns = 4, lines = TRUE, points = FALSE),
        xlab = "Effort to arise from stool (Borg scale)")
print(fm2 <- lmer(effort ~ Type + (1|Subject), ergoStool), corr = FALSE)
print(fm2a <- lmer(effort ~ 0 + Type + (1|Subject), ergoStool), corr = FALSE)
(fm3 <- lmer(effort ~ 1 + (1| Type) + (1|Subject), ergoStool))
invisible(rr <- dotplot(ranef(fm3, post = TRUE)))
rr$Type
rr$Subject

str(sch1 <- unique(subset(MathAchieve, select = c(School, MEANSES))))
xtabs(~ xtabs(~ School, MathAchieve))
densityplot(~ xtabs(~ School, MathAchieve))
(fm4 <- lmer(MathAch ~ Sex + Minority + MEANSES + (1|School), MathAchieve))
dotplot(ranef(fm4, post = TRUE), strip = FALSE)
qqmath(ranef(fm4, post = TRUE), strip = FALSE)
(fm5 <- lmer(MathAch ~ Sex * Minority + MEANSES + (1|School), MathAchieve))

xyplot(weight ~ Lsize, RatPupWeight, groups = sex, type = c("g","p","smooth"))
xyplot(weight ~ log(Lsize), RatPupWeight, groups = sex, type = c("g","p","smooth"))
xyplot(weight ~ Lsize, RatPupWeight, groups = sex,
       scales = list(x = list(log = 2)), type = c("g","p","smooth"),
       auto.key = list(columns = 2, lines = TRUE, points = FALSE))
xyplot(weight ~ Lsize, RatPupWeight, groups = Treatment,
       scales = list(x = list(log = 2)), type = c("g","p","smooth"),
       auto.key = list(columns = 3, lines = TRUE, points = FALSE))
xtabs(~ Treatment + Litter, RatPupWeight, sparse = TRUE)
(fm6 <- lmer(weight ~ sex + Lsize + Treatment + (1|Litter), RatPupWeight))
