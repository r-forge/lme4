###################################################
### chunk number 1: preliminaries
###################################################
#line 19 "/home/bates/Documents/slides/2011-01-11-Madison/4Theory.Rnw"
options(width=69,show.signif.stars=FALSE)
library(lattice)
lattice.options(default.theme = function() standard.theme())
#lattice.options(default.theme = function() standard.theme(color=FALSE))
library(Matrix)
library(lme4a)
if (file.exists("classroom.rda")) {
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

simplemer <- function(flist, y, X, REML = TRUE, super = FALSE)
{
    ## check arguments for consistency
    stopifnot(is.numeric(y),
              is.matrix(X) | is(X, "Matrix"),
              (n <- length(y)) > 0,
              nrow(X) == n,
              is.list(flist),
              length(flist) > 0,
              all(sapply(flist, is.factor)),
              all(sapply(flist, length) == n))
    super <- as.logical(super)[1]

    rho <- new.env(parent = emptyenv()) # create an empty environment
    rho$y <- y                    # store arguments and derived values
    rho$X <- X
    chol(rho$XtX <- crossprod(X))       # check for full column rank
    rho$REML <- as.logical(REML)[1]

    rho$Zt <- do.call(rBind, lapply(flist, as, "sparseMatrix"))
    rho$nlev <- sapply(flist, function(x) length(levels(factor(x))))
    rho$L <- Cholesky(tcrossprod(rho$Zt), LDL = FALSE, Imult = 1, super = super)
    return(rho)
}


###################################################
### chunk number 2: PenicillinL
###################################################
#line 284 "/home/bates/Documents/slides/2011-01-11-Madison/4Theory.Rnw"
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
#line 297 "/home/bates/Documents/slides/2011-01-11-Madison/4Theory.Rnw"
print(image(tcrossprod(Ut), xlab = NULL, ylab = NULL,
            sub = expression(Lambda*"'Z"*"'Z"*Lambda+"I")),
      pos = c(0,0,0.47,1), more = TRUE)
print(image(L, xlab = NULL, ylab = NULL, sub = "L"),
      pos = c(0.47,0,1,1))


###################################################
### chunk number 4: revChol
###################################################
#line 322 "/home/bates/Documents/slides/2011-01-11-Madison/4Theory.Rnw"
Zt <- do.call(rBind, lapply(flist[2:1], as, "sparseMatrix"))
Lambda <- Diagonal(x = rep.int(theta[2:1], nlev[2:1]))
Ut <- crossprod(Lambda, Zt)
Lnoperm <- Cholesky(tcrossprod(Ut), perm = FALSE, LDL = FALSE, Imult = 1)
Lperm <- Cholesky(tcrossprod(Ut), LDL = FALSE, Imult = 1)
sapply(lapply(list(L, Lnoperm, Lperm), as, "sparseMatrix"), nnzero)


###################################################
### chunk number 5: Reversedfactorimages
###################################################
#line 334 "/home/bates/Documents/slides/2011-01-11-Madison/4Theory.Rnw"
print(image(Lnoperm, xlab = NULL, ylab = NULL, sub = "Lnoperm"),
      split = c(1,1,2,1), more = TRUE)
print(image(Lperm, xlab = NULL, ylab = NULL, sub = "Lperm"),
      split = c(2,1,2,1))


###################################################
### chunk number 6: classroomL
###################################################
#line 361 "/home/bates/Documents/slides/2011-01-11-Madison/4Theory.Rnw"
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
#line 376 "/home/bates/Documents/slides/2011-01-11-Madison/4Theory.Rnw"
print(image(tcrossprod(Ut), xlab = NULL, ylab = NULL, 
            sub = expression(Lambda*"'Z"*"'Z"*Lambda+"I")),
      split = c(1,1,2,1), more = TRUE)
print(image(L, xlab = NULL, ylab = NULL, sub = "L"),
      split = c(2,1,2,1))


###################################################
### chunk number 8: fm1
###################################################
#line 602 "/home/bates/Documents/slides/2011-01-11-Madison/4Theory.Rnw"
invisible(lmer(mathgain ~ mathkind + minority + ses + (1|classid) + (1|schoolid), classroom, verbose = 1, REML = FALSE))


