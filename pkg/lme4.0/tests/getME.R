library(lme4.0)
#### tests of getME()

###  are names correct? --------------
if(getRversion() < "2.15")
    paste0 <- function(...) paste(..., sep = '')
hasInms <- function(x) grepl("(Intercept", names(x), fixed=TRUE)
matchNms <- function(fm, PAR) {
    stopifnot(is.character(vnms <- names(fm@flist)))
    mapply(grepl, paste0("^", vnms), names(PAR))
}
chkIMod <- function(fm) {## check "intercept only" model
    b1 <- getME(fm,"beta")
    f1 <- fixef(fm)
    stopifnot(hasInms(f1), f1 == b1,
              hasInms(t1 <- getME(fm,"theta")), matchNms(fm, t1))
}

fm1 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
chkIMod(fm1)

fm2 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake)
stopifnot(fixef(fm2) == getME(fm2,"beta"))
getME(fm2,"theta")

getME(fm3 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy),
      "theta")
getME(fm4 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy),
      "theta")

## internal consistency check ensuring that all allowed 'name's work (and are not empty):
nmME <- eval(formals(getME)$name)
## lme4.0 does not yet support all:  remove current exceptions:
(nmME <- nmME[!(nmME %in% c("u", "REML"))])

chkMEs <- function(fm, nms) {
    stopifnot(is.character(nms))
    str(parts <- sapply(nms, getME, object = fm, simplify=FALSE))
    isN <- sapply(parts, is.null)
    stopifnot(identical(names(isN), nms), !any(isN))
}

chkMEs(fm1, nmME)
chkMEs(fm2, nmME)
chkMEs(fm3, nmME)
chkMEs(fm4, nmME)

isREML(fm1)

L <- as(fm1@L,"Matrix")
Z <- getME(fm1,"Z")
A <- fm1@A
dim(L)
dim(Z)
dim(A)
