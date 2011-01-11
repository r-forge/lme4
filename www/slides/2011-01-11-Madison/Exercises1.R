###################################################
### chunk number 1: initial
###################################################
#line 20 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises1.Rnw"
options(width=80, show.signif.stars = FALSE)
library(lattice)
lattice.options(default.theme = standard.theme(color = FALSE))
set.seed(123454321)


###################################################
### chunk number 2: 
###################################################
#line 42 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises1.Rnw"
classroom <- read.csv("http://www-personal.umich.edu/~bwest/classroom.csv")


###################################################
### chunk number 3: strclass1
###################################################
#line 46 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises1.Rnw"
foo <- capture.output(str(classroom))


###################################################
### chunk number 4: strclass eval=FALSE
###################################################
## #line 53 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises1.Rnw"
## str(classroom)


###################################################
### chunk number 5: strclass2
###################################################
#line 56 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises1.Rnw"
cat(paste(foo[1:5],collapse="\n"),"\n")


###################################################
### chunk number 6: sexchange
###################################################
#line 59 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises1.Rnw"
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
#line 80 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises1.Rnw"
summary(classroom$sex)


###################################################
### chunk number 8: subsetsummary
###################################################
#line 90 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises1.Rnw"
summary(subset(classroom, select = c(sex, minority, childid, classid, schoolid)))


###################################################
### chunk number 9: saveload
###################################################
#line 97 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises1.Rnw"
save(classroom, file = "classroom.rda")
rm(classroom)
exists("classroom")
load("classroom.rda")


###################################################
### chunk number 10: checkstr
###################################################
#line 103 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises1.Rnw"
str(classroom)


###################################################
### chunk number 11: checkchildid
###################################################
#line 112 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises1.Rnw"
all(row.names(classroom) == classroom$childid)


###################################################
### chunk number 12: childid
###################################################
#line 118 "/home/bates/Documents/slides/2011-01-11-Madison/Exercises1.Rnw"
classroom$childid <- NULL
names(classroom)


