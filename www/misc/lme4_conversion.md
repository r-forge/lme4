# lme4 release guide

## Version numbering

As previously announced on the mailing list, we will shortly be releasing a new version of `lme4`, a descendant of the previous development version `lme4Eigen`. For users who do not access any internal structures, there will be few backward-incompatible changes.

* the version of `lme4` currently on R-forge (currently version 0.999902344-0, to be released as version 1.0 (!)) should be used for new projects
* The current CRAN version (0.999375-42) will be replaced by a nearly identical version called `lme4.0` (currently version 0.9999-2; the only backward-incompatible change in this version is a fix to the AGQ code contributed by Wayne Zhang).  `lme4.0` is a maintenance version and will only be changed to fix documented bugs.
* all other versions (`lme4a`, `lme4b`, `lme4Eigen` from R-forge) are deprecated.

## Changes in behavior
* Because the internal computational machinery has changed, results from the newest version of `lme4` will not be numerically identical to those from previous versions.  For reasonably well-defined fits, they will be extremely close (within numerical tolerances of 1e-4 or so), but for unstable or poorly-defined fits the results may change, and very unstable fits may fail when they (apparently) succeeded with previous versions. Similarly, some fits may be slower with the new version, although on average the new version should be faster and more stable. There are more numerical tuning options available than before (see below); non-default settings may restore the speed and/or ability to fit a particular model without an error.
* In the past, `lmer` automatically called `glmer` when `family` was specified. It still does so, but now warns the user that they should preferably use `glmer` directly
* `VarCorr` returns its results in the same format as before (as a list of variance-covariance matrices with `correlation` and `stddev` attributes, plus a `sc` attribute giving the residual standard deviation/scale parameter when appropriate), but prints them in a different (nicer) way

## Other user-visible changes
* methods such as `fixef` no longer conflict with the `nlme` package
* `[gn]lmer` now produces objects of class `merMod` rather than class `mer` as before
* `[gn]lmer` now has an `optimizer` argument. `"Nelder_Mead"` is the default for `[n]lmer`, while a combination of `"bobyqa"` (an alternative derivative-free method) and `"Nelder_Mead"` is the default for `glmer`; to use the `nlminb` optimizer as in the old version of lme4, you can use `optimizer="optimx"` with `control=list(method="nlminb")` (you will need the `optimx` package to be installed and loaded). See the help pages for details.
* (FIXME: describe residuals)

## New features
* a general-purpose `getME()` accessor method has been used to allow extraction of a wide variety of components of a mixed-model fit; this has been backported to `lme4.0` for compatibility
* `bootMer`, a framework for obtaining parameter confidence intervals by parametric bootstrapping
* `plot` methods similar to those from the `nlme` package (although without `augPred`)
* a `predict` method, allowing a choice of which random effects are included in the prediction
* profiling (and profile confidence intervals) for `lmer` and `glmer` results
* `nAGQ=0`, an option to do fast (but inaccurate) fitting of GLMMs

## Still non-existent features
* Automatic MCMC sampling based on the fit turns out to be very difficult to implement in a way that is really broadly reliable and robust; `mcmcsamp` will not be implemented in the near future
* "R-side" structures (within-block correlation and heteroscedasticity) are not on the current timetable

## For package writers:

[Current package compatibilty test results][pkgtest]

[pkgtest]: ./lme4_compat_report.html

* `lme4`-old and `lme4.0` produces `mer` objects, `lme4`-new produces `merMod` objects
* you can distinguish `lme4`-old from `lme4`-new via package version; the last old-style version of `lme4` on CRAN is 0.999375-42, so anything after that is `lme4`-new (the current version on <http://lme4.r-forge.r-project.org/repos> is 0.999902344-0)
* so you can test e.g. if `packageVersion("lme4")<="0.999375-43"` (yes, you do want the quotation marks; package versions are weird objects in R)
* `getME(.,.)` should often get the components you want
* `lme4` now uses S3 rather than S4 methods in many places. That makes some things easier ...

### Things that won't work:
* direct extraction of slots via `@` (use `getME()` instead)
* methods that depend on `lme4` producing objects of class `mer` (write new methods for class `merMod`; `lme4` now has `isLMM()`, `isGLMM()`, `isNLMM()` that should help you distinguish different types of model if you need to)
* the `method` argument is no longer used in `glmer` (`nAGQ=1` vs `nAGQ>1` specifies whether to use Laplace or 
* `expandSlash` no longer exists (*although it does exist within the `findbars` function: could be reconstituted??*)
* because S4 methods are used less, and (S4) *reference* classes are used, considerably fewer methods are exported, but they are generally available as reference class method "slots"

