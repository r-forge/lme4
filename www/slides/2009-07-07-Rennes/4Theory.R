###################################################
### chunk number 1: preliminaries
###################################################
options(width=69,show.signif.stars=FALSE)
library(lattice)
lattice.options(default.theme = function() standard.theme())
lattice.options(default.theme = function() standard.theme(color=FALSE))
library(Matrix)
library(lme4)
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

profDev <- function(rho, theta) {
    stopifnot(is.numeric(theta), length(theta) == length(rho$nlev))
    Ut <- crossprod(Diagonal(x = rep.int(theta, rho$nlev)), rho$Zt)
    L <- update(rho$L, Ut, mult = 1)
    cu <- solve(L, solve(L, Ut %*% rho$y, sys = "P"), sys = "L")
    RZX <- solve(L, solve(L, Ut %*% rho$X, sys = "P"), sys = "L")
    RX <- chol(rho$XtX - crossprod(RZX))
    cb <- solve(t(RX), crossprod(rho$X, rho$y) - crossprod(RZX, cu))
    beta <- solve(RX, cb)
    u <- solve(L, solve(L, cu - RZX %*% beta, sys = "Lt"), sys = "Pt")
    fitted <- as.vector(crossprod(Ut, u) + rho$X %*% beta)
    prss <- sum(c(rho$y - fitted, as.vector(u))^2) # penalized residual sum of squares
    n <- length(fitted);  p <- ncol(RX)
    if (rho$REML) return(determinant(L)$mod + 2*determinant(RX)$mod +
                         (n-p)*(1+log(2*pi*prss/(n-p))))
    determinant(L)$mod + n * (1 + log(2*pi*prss/n))}


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


###################################################
### chunk number 8: fm1
###################################################
invisible(lmer(mathgain ~ mathkind + minority + ses + (1|classid) + (1|schoolid), classroom, verbose = 1, REML = FALSE))


###################################################
### chunk number 9: profDev1
###################################################
ls.str(rho <- with(classroom, simplemer(list(classid, schoolid), mathgain, model.matrix(~ mathkind + minority + ses), REML = FALSE)))
invisible(nlminb(c(0.836158, 0.489669), function(x) profDev(rho, x), lower = c(0,0), control = list(trace = 1)))


