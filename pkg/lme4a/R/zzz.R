.NameSpace <- environment()

.onLoad <- function(libname, pkgname)
{
     mod <- Module("lme4")
     populate(mod, .NameSpace)
}
