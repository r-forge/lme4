# lme4 updating guide

## For package writers:

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

## For package writers and regular people:
* accessor functions and methods should all work the same (slight differences in printing format for `VarCorr`)
* there is now an `optimizer` argument. `"Nelder_Mead"` is the default, other options are `"bobyqa"` (an alternative derivative-free method); to use the `nlminb` optimizer as in the old version of lme4, you can use `optimizer="optimx"` with `control=list(method="nlminb")` (you will need the `optimx` package to be installed and loaded)
* (residuals)

## New features
* `bootMer`
* profiling
* `nAGQ=0`
* choice of optimizer, including user-specified

## User-visible changes
* `glmer` not automatically called when `family` is specified in `lmer`?
* methods such as `fixef` no longer conflict between `nlme` and `lme4`
