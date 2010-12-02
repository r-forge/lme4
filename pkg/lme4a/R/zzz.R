.NameSpace <- environment()

.onLoad <- function(libname, pkgname)
{
    ## For Matrix API change (Oct.2009) - silence the warning:
    assign("det_CHMfactor.warn", FALSE, envir = Matrix:::.MatrixEnv)
    mod <- Module("lme4")
    populate(mod, .NameSpace)
}
