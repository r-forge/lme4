###################################################
### chunk number 1: preliminaries
###################################################
options(width=69,show.signif.stars=FALSE,str=strOptions(strict.width = "cut"))


###################################################
### chunk number 2: readcsv eval=FALSE
###################################################
## mydata <- read.csv(file.choose())


###################################################
### chunk number 3: readcsv eval=FALSE
###################################################
## mydata <- read.delim(file.choose())


###################################################
### chunk number 4: Formaldehyde
###################################################
str(Formaldehyde)
summary(Formaldehyde)
Formaldehyde


###################################################
### chunk number 5: InsectSprays
###################################################
str(InsectSprays)
summary(InsectSprays)
head(InsectSprays)


###################################################
### chunk number 6: saverestore
###################################################
sprays <- InsectSprays
save(sprays, file = "sprays.rda")
rm(sprays)
ls.str()
load("sprays.rda")
names(sprays)


###################################################
### chunk number 7: dollarop
###################################################
Formaldehyde$carb


###################################################
### chunk number 8: dollaropleft
###################################################
sprays$sqrtcount <- sqrt(sprays$count)
names(sprays)


###################################################
### chunk number 9: dollaropleftNULL
###################################################
sprays$sqrtcount <- NULL
names(sprays)


###################################################
### chunk number 10: formalfoo
###################################################
Formaldehyde$carb * Formaldehyde$optden
with(Formaldehyde, carb * optden)


###################################################
### chunk number 11: within
###################################################
sprays <- within(sprays, sqrtcount <- sqrt(count))
str(sprays)


###################################################
### chunk number 12: sprays
###################################################
str(sprays <- within(InsectSprays, spray <- as.integer(spray)))
str(sprays <- within(sprays, spray <- factor(spray, labels = LETTERS[1:6])))


###################################################
### chunk number 13: sprayA
###################################################
str(sprayA <- subset(sprays, spray == "A"))


###################################################
### chunk number 14: xtabssprays
###################################################
xtabs( ~ spray, sprayA)
xtabs( ~ spray, sprayA, drop = TRUE)


###################################################
### chunk number 15: sprayAdrop
###################################################
str(sprayA <- within(sprayA, spray <- factor(spray)))
xtabs( ~ spray, sprayA)


###################################################
### chunk number 16: sprayDEF
###################################################
str(sprayDEF <- subset(sprays, spray %in% c("D","E","F")))


###################################################
### chunk number 17: unstack
###################################################
str(unstack(InsectSprays))


###################################################
### chunk number 18: <classprep
###################################################
if (file.exists("classroom.rda")) {
  load("classroom.rda")
} else {
    classroom <- within(read.csv("http://www-personal.umich.edu/~bwest/classroom.csv"),
                        schoolid <- factor(schoolid))
}


###################################################
### chunk number 19: clasuniq
###################################################
str(unique(subset(classroom, select = c(schoolid,housepov))))


###################################################
### chunk number 20: <cleanup
###################################################
rm(classroom, sprayA, sprayDEF, sprays)


