\name{Dyestuff}
\alias{Dyestuff}
\docType{data}
\title{Yield of dyestuff by batch}
\description{
  Yield of dyestuff (Naphthalene Black 12B) from 5 different
  preparations from each of 6 different batchs of an intermediate
  product (H-acid).
}
\usage{data(Dyestuff)}
\format{
  A data frame with 30 observations on the following 2 variables.
  \describe{
    \item{\code{Batch}}{a factor indicating the batch of the
      intermediate product from which the preparation was created.}
    \item{\code{Yield}}{the yield of dyestuff from the preparation
      (grams of standard color).}
  }
}
\details{
 The data are described in Davies and Goldsmith (1972) as coming from
 "an investigation to find out how much the variation from batch to
 batch in the quality of an intermediate product (H-acid) contributes to
 the variation in the yield of the dyestuff (Naphthalene Black 12B) made
 from it.  In the experiment six samples of the intermediate,
 representing different batches of works manufacture, were obtained, and
 five preparations of the dyestuff were made in the laboratory from each
 sample.  The equivalent yield of each preparation as grams of standard
 colour was determined by dye-trial".
}
\source{
 O.L. Davies and P.L. Goldsmith (eds), \emph{Statistical Methods in
   Research and Production, 4th ed.}, Oliver and Boyd, (1972), section
   6.4
}
%\references{}
\examples{
str(Dyestuff)
dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff,
        ylab = "Batch", jitter.y = TRUE, aspect = 0.3,
        type = c("p", "a"))
(fm1 <- lmer(Yield ~ 1|Batch, Dyestuff))
}
\keyword{datasets}
