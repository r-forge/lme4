\name{Multilocation}
\alias{Multilocation}
\title{A multilocation trial}
\description{
  The \code{Multilocation} data frame has 108 rows and 7 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{obs}{
      a numeric vector
    }
    \item{Location}{
      an ordered factor with levels
      \code{B} < \code{D} < \code{E} < \code{I} < \code{G} <
      \code{A} < \code{C} < \code{F} < \code{H}
    }
    \item{Block}{
      a factor with levels \code{1} to \code{3} 
    }
    \item{Trt}{
      a factor with levels \code{1} to \code{4}
    }
    \item{Adj}{
      a numeric vector
    }
    \item{Fe}{
      a numeric vector
    }
    \item{Grp}{
      an ordered factor with levels
      \code{B/1} < \code{B/2} < \code{B/3} < \code{D/1} <
      \code{D/2} < \code{D/3} < \code{E/1} < \code{E/2} <
      \code{E/3} < \code{I/1} < \code{I/2} < \code{I/3} <
      \code{G/1} < \code{G/2} < \code{G/3} < \code{A/1} <
      \code{A/2} < \code{A/3} < \code{C/1} < \code{C/2} <
      \code{C/3} < \code{F/1} < \code{F/2} < \code{F/3} <
      \code{H/1} < \code{H/2} < \code{H/3}
    }
  }
}
\source{
  Littel, R. C., Milliken, G. A., Stroup, W. W., and Wolfinger,
  R. D. (1996), \emph{SAS System for Mixed Models}, SAS Institute
  (Data Set 2.8.1).
}
\examples{
str(Multilocation)
if (require("lme4", quietly = TRUE, character = TRUE)) {
  options(contrasts = c(unordered = "contr.SAS", ordered = "contr.poly"))
  ### Create a Block %in% Location factor
  Multilocation$Grp <- with(Multilocation, Block:Location)
  print(fm1Mult <- lmer(Adj ~ Location * Trt + (1|Grp), Multilocation))
  print(anova(fm1Mult))
  print(fm2Mult <- lmer(Adj ~ Location + Trt + (1|Grp), Multilocation))
  print(fm3Mult <- lmer(Adj ~ Location + (1|Grp), Multilocation))
  print(fm4Mult <- lmer(Adj ~ Trt + (1|Grp), Multilocation))
  print(fm5Mult <- lmer(Adj ~ 1 + (1|Grp), Multilocation))
  print(anova(fm2Mult))
  print(anova(fm1Mult, fm2Mult, fm3Mult, fm4Mult, fm5Mult))
  ### Treating the location as a random effect
  print(fm1MultR <- lmer(Adj ~ Trt + (1|Location/Trt) + (1|Grp), Multilocation))
  print(anova(fm1MultR))
  print(fm2MultR <- lmer(Adj ~ Trt + (Trt - 1|Location) + (1|Block),
                         Multilocation))
  print(anova(fm1MultR, fm2MultR))
\dontrun{
intervals(fm1MultR)
}
}
}
\keyword{datasets}
