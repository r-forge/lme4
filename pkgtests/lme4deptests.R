## code for testing lme4 downstream packages
##
## include all downstream packages from CRAN, r-forge:
## packages to check, loaded from package-dependency scan

source("http://developer.r-project.org/CRAN/Scripts/depends.R")
rr <- reverse_dependencies_with_maintainers("lme4")
source("lme4depfuns.R")  ## component/utility functions
## packages to skip (because dependencies don't install etc.
## skippkgs <- "polytomous"
skippkgs <- character(0)

## * must export R_LIBS_SITE=./library before running R CMD BATCH
##   and  make sure that .R/check.Renviron is set
##   (this is done in the 'runtests' script)

##  (repeated from lme4depfuns.R):
##   FIXME: check for/document local version of tarball more recent than R-forge/CRAN versions
##  currently just tries to find most recent version and check it, but could check all versions?
##  FIXME: consistent implementation of checkdir

## FIXME: set up an appropriate makefile structure for this ? (a little tricky if it also depends on
##   checking CRAN/R-forge versions?
##  might to be able to use update.packages() ...
  
## FIXME (?)/warning: need to make sure that latest/appropriate version of lme4 is installed locally ...

## FIXME: why are R2admb, RLRsim, sdtalt, Zelig not getting checked?

testdir <- getwd()

## directory for tarballs to check
tarballdir <- file.path(testdir,"tarballs")
## directory for ?? (at least current lme4 version ...)
libdir <-     file.path(testdir,"library")

suppressWarnings(rm(list=c("availCRAN","availRforge"))) ## clean up

## require(tools)

## make directories ...
dir.create(tarballdir,showWarnings=FALSE)
dir.create(libdir,showWarnings=FALSE)
## want to install additional dependencies etc. out of the way
## to keep original installed base clean, but this may not be feasible
## it would be nice to use tools:::testInstalledPackages(), but I may simply
##  have to do R CMD check

pkgnames <- rr[,"Package"]
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

