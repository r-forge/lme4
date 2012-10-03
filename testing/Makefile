#####
Cbpp.pdf: Cbpp_lme4_bobyqa.RData Cbpp_lme4_NM.RData Cbpp_lme4.0.RData Cbpp_glmmML.RData Cbpp_glmmADMB.RData Cbpp.Rnw
	Rscript -e "library(knitr); knit2pdf(\"Cbpp.Rnw\")"

Generic_lme4.R: miscfuns.R

Cbpp_lme4_bobyqa.RData: Cbpp_lme4_bobyqa.R Generic_lme4.R
	R CMD BATCH --vanilla Cbpp_lme4_bobyqa.R

Cbpp_lme4_NM.RData: Cbpp_lme4_NM.R Generic_lme4.R
	R CMD BATCH --vanilla Cbpp_lme4_NM.R

Cbpp_lme4.0.RData: Cbpp_lme4.0.R Generic_lme4.0.R
	R CMD BATCH --vanilla Cbpp_lme4.0.R

Cbpp_glmmML.RData: Cbpp_glmmML.R
	R CMD BATCH --vanilla Cbpp_glmmML.R

Cbpp_glmmADMB.RData: Cbpp_glmmadmb.R
	R CMD BATCH --vanilla Cbpp_glmmadmb.R

###############
Bacteria.pdf: Bacteria_lme4_bobyqa.RData Bacteria_lme4_NM.RData Bacteria_lme4.0.RData Bacteria_glmmML.RData Bacteria_glmmADMB.RData Bacteria.Rnw
	Rscript -e "library(knitr); knit2pdf(\"Bacteria.Rnw\")"

Bacteria_lme4_bobyqa.RData: Bacteria_lme4_bobyqa.R Generic_lme4.R
	R CMD BATCH --vanilla Bacteria_lme4_bobyqa.R

Bacteria_lme4_NM.RData: Bacteria_lme4_NM.R Generic_lme4.R
	R CMD BATCH --vanilla Bacteria_lme4_NM.R

Bacteria_lme4.0.RData: Bacteria_lme4.0.R Generic_lme4.0.R
	R CMD BATCH --vanilla Bacteria_lme4.0.R

Bacteria_glmmML.RData: Bacteria_glmmML.R
	R CMD BATCH --vanilla Bacteria_glmmML.R

Bacteria_glmmADMB.RData: Bacteria_glmmadmb.R
	R CMD BATCH --vanilla Bacteria_glmmadmb.R

################
Contraception_lme4e_bobyqa.RData: Contraception_lme4e_bobyqa.R
	R CMD BATCH --vanilla Contraception_lme4_bobyqa.R

Contraception_lme4e_NM.RData: Contraception_lme4e_NM.R
	R CMD BATCH --vanilla Contraception_lme4e_NM.R

