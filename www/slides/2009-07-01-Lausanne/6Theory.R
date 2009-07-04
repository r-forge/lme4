###################################################
### chunk number 1: preliminaries
###################################################
options(width=69,show.signif.stars=FALSE)
library(lattice)
lattice.options(default.theme = function() standard.theme())
#lattice.options(default.theme = function() standard.theme(color=FALSE))
library(Matrix)
#library(lme4)
data(Dyestuff, package = "lme4")
data(Penicillin, package = "lme4")
data(Pastes, package = "lme4")
if (exists("classroom.rda")) {
    load("classroom.rda")
} else {
    classroom <- within(read.csv("http://www-personal.umich.edu/~bwest/classroom.csv"),
                    {
                        classid <- factor(classid)
                        schoolid <- factor(schoolid)
                        sex <- factor(sex, labels = c("M","F"))
                        minority <- factor(minority, labels = c("N", "Y"))
                    })
    save(classroom, file = "classroom.rda")
}


###################################################
### chunk number 2: PenicillinL
###################################################
flist <- subset(Penicillin, select = c(plate, sample))
Zt <- do.call(rBind, lapply(flist, as, "sparseMatrix"))
(nlev <- sapply(flist, function(f) length(levels(factor(f)))))
theta <- c(1.2, 2.1)
Lambda <- Diagonal(x = rep.int(theta, nlev))
Ut <- crossprod(Lambda, Zt)
str(L <- Cholesky(tcrossprod(Ut), LDL = FALSE, Imult = 1))


###################################################
### chunk number 3: Penicillinimage
###################################################
print(image(tcrossprod(Ut), xlab = NULL, ylab = NULL, sub = "U'U+I"),
      pos = c(0,0,0.47,1), more = TRUE)
print(image(L, xlab = NULL, ylab = NULL, sub = "L"),
      pos = c(0.47,0,1,1))


###################################################
### chunk number 4: revChol
###################################################
Zt <- do.call(rBind, lapply(flist[2:1], as, "sparseMatrix"))
Lambda <- Diagonal(x = rep.int(theta[2:1], nlev[2:1]))
Ut <- crossprod(Lambda, Zt)
Lnoperm <- Cholesky(tcrossprod(Ut), perm = FALSE, LDL = FALSE, Imult = 1)
Lperm <- Cholesky(tcrossprod(Ut), LDL = FALSE, Imult = 1)
sapply(lapply(list(L, Lnoperm, Lperm), as, "sparseMatrix"), nnzero)


###################################################
### chunk number 5: Reversedfactorimages
###################################################
print(image(Lnoperm, xlab = NULL, ylab = NULL, sub = "Lnoperm"),
      split = c(1,1,2,1), more = TRUE)
print(image(Lperm, xlab = NULL, ylab = NULL, sub = "Lperm"),
      split = c(2,1,2,1))


###################################################
### chunk number 6: classroomL
###################################################
Zt <- do.call(rBind, lapply(flist <- subset(Pastes,,c(sample, batch)),
                            as, "sparseMatrix"))
nlev <- sapply(flist, function(f) length(levels(factor(f))))
theta <- c(0.4, 0.5)
Lambda <- Diagonal(x = rep.int(theta, nlev))
Ut <- crossprod(Lambda, Zt)
L <- Cholesky(tcrossprod(Ut), LDL = FALSE, Imult = 1)
str(L@perm)


###################################################
### chunk number 7: Pastesimage
###################################################
print(image(tcrossprod(Ut), xlab = NULL, ylab = NULL, sub = "U'U+I"),
      split = c(1,1,2,1), more = TRUE)
print(image(L, xlab = NULL, ylab = NULL, sub = "L"),
      split = c(2,1,2,1))


