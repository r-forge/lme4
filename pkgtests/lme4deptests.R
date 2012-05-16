## code for testing lme4 downstream packages
##
## include all downstream packages from CRAN, r-forge:
## packages to check, loaded from package-dependency scan
X <- read.csv("lme4_depends.csv",as.is=TRUE,fileEncoding="UTF8")
source("lme4depfuns.R")
skippkgs <- "polytomous"

## * must export R_LIBS_SITE=./library before running R CMD BATCH?  (why?)
## * make sure that .R/check.Renviron is set, e.g.
##   mv /home/bolker/.R/check.Renviron.old  ~/.R/check.Renviron
##     ~/.R/check.Renviron:
##        R_LIBS_SITE=/home/bolker/R/pkgs/lme4/testpkgs/library

##  (repeated from lme4depfuns.R):
##    FIXME: check for/document local version of tarball more recent than R-forge/CRAN versions
##    FIXME: consistent implementation of checkdir

## FIXME: set up an appropriate makefile structure for this ? (a little tricky if it also depends on
##   checking CRAN/R-forge versions?

## FIXME (?)/warning: need to make sure that latest version of lme4 is installed in libdir ...

## FIXME: why are R2admb, RLRsim, sdtalt, Zelig not getting checked?

testdir <- "~/R/pkgs/lme4/pkgtests"
tarballdir <- file.path(testdir,"tarballs")
libdir <-     file.path(testdir,"library")

suppressWarnings(rm(list=c("availCRAN","availRforge"))) ## clean up
options(repos="http://probability.ca/cran")

require(tools)
## assume wd is test dir, already
if (FALSE) {
    setwd(testdir)
}

## auto-generated list -- cleaned up by hand (removed 'lme4', fixed R2STATS maintainer syntax,
##   fixed comma in LMERconveniencefunctions author)

dir.create(tarballdir,showWarnings=FALSE)
dir.create(libdir,showWarnings=FALSE)
## want to install additional dependencies etc. out of the way
## to keep original installed base clean, but this may not be feasible
## it would be nice to use tools:::testInstalledPackages(), but I may simply
##  have to do R CMD check

pkgnames <- X$X
names(pkgnames) <- pkgnames ## so results are named
pkgnames <- pkgnames[!pkgnames %in% skippkgs]
require(parallel)
## FIXME (maybe): mclapply doesn't work on Windows ?
testresults <- mclapply(pkgnames,function(x) try(checkPkg(x,verbose=TRUE)))
skipresults <- mclapply(skippkgs,function(x) try(checkPkg(x,skip=TRUE,verbose=TRUE)))
testresults <- c(testresults,skipresults)

save("testresults",file="lme4tests_out.RData")

if (FALSE) {
    ## playing with results
    load("lme4tests_out.RData")
    checkPkg("HSAUR2")
}
## names(testresults) <- X$X  ## should be obsolete after next run
genReport(X,testresults)

