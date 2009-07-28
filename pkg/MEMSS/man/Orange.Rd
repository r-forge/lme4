\name{Orange}
\docType{data}
\alias{Orange}
\title{Growth of orange trees}
\description{
  The \code{Orange} data frame has 35 rows and 3 columns of records of
  the growth of orange trees.
}
\usage{Orange}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Tree}{
      a factor indicating the tree on which the measurement is
      made.
    }
    \item{age}{
      a numeric vector giving the age of the tree (days since 1968/12/31)
    }
    \item{circumference}{
      a numeric vector of trunk circumferences (mm).  This is probably
      \dQuote{circumference at breast height}, a standard measurement in
      forestry.
    }
  }
}
\source{
  Draper, N. R. and Smith, H. (1998), \emph{Applied Regression Analysis
    (3rd ed)}, Wiley (exercise 24.N).

  Pinheiro, J. C. and Bates, D. M. (2000) \emph{Mixed-effects Models
    in S and S-PLUS}, Springer.
}
\examples{
xyplot(circumference ~ age, Orange, groups = Tree, type = c("g", "b"),
       auto.key = list(space = "right", lines = TRUE), aspect = "xy",
       xlab = "Age (days since 1968/12/31)", ylab = "Circumference (mm)")
\dontrun{
m1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
            Orange, verbose = TRUE,
            start = c(Asym = 190, xmid = 730, scal = 350))
.Call("mer_optimize", m1, 1L, 1L, PACKAGE = "lme4")
print(m1)
ranef(m1)
}
}
\keyword{datasets}
