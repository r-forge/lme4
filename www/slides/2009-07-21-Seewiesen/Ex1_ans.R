classroom <- read.csv("http://www-personal.umich.edu/~bwest/classroom.csv")
str(classroom)
summary(classroom$sex)

# after some more experimentation

classroom <- within(read.csv("http://www-personal.umich.edu/~bwest/classroom.csv"),
                {
                    sex <- factor(sex, labels = c("M","F"))
                    minority <- factor(minority, labels = c("N","Y"))
                    classid <- factor(classid)
                    schoolid <- factor(schoolid)
                    childid <- factor(childid)
                })

summary(subset(classroom, select = c(sex, minority, childid, classid, schoolid)))
save(classroom, file = "classroom.rda")
all(row.names(classroom) == classroom$childid)
classroom$childid <- NULL
names(classroom)

str(school <- unique(subset(classroom, select = c(schoolid, housepov))))
str(class <- unique(subset(classroom, select = c(classid, schoolid,
                                      yearstea, mathknow, mathprep))))
str(classes <- merge(class, school))
