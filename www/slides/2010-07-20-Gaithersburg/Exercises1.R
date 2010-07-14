###################################################
### chunk number 1: initial
###################################################
options(width=76, show.signif.stars = FALSE)
library(lattice)
lattice.options(default.theme = standard.theme(color = FALSE))
set.seed(123454321)


###################################################
### chunk number 2: 
###################################################
classroom <- read.csv("http://www-personal.umich.edu/~bwest/classroom.csv")


###################################################
### chunk number 3: strclass1
###################################################
foo <- capture.output(str(classroom))


###################################################
### chunk number 4: strclass eval=FALSE
###################################################
## str(classroom)


###################################################
### chunk number 5: strclass2
###################################################
cat(paste(foo[1:5],collapse="\n"),"\n")


###################################################
### chunk number 6: sexchange
###################################################
classroom <- within(classroom,
                {
                    sex <- factor(sex, labels = c("M","F"))
                    minority <- factor(minority, labels = c("N","Y"))
                    classid <- factor(classid)
                    schoolid <- factor(schoolid)
                    childid <- factor(childid)
                })


###################################################
### chunk number 7: sexsummary
###################################################
summary(classroom$sex)


###################################################
### chunk number 8: subsetsummary
###################################################
summary(subset(classroom, select = c(sex, minority, childid, classid, schoolid)))


###################################################
### chunk number 9: saveload
###################################################
save(classroom, file = "classroom.rda")
rm(classroom)
exists("classroom")
load("classroom.rda")


###################################################
### chunk number 10: checkstr
###################################################
str(classroom)


###################################################
### chunk number 11: checkchildid
###################################################
all(row.names(classroom) == classroom$childid)


###################################################
### chunk number 12: childid
###################################################
classroom$childid <- NULL
names(classroom)


