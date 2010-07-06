setAs("ddenseModelMatrix", "predModule",
      function(from)
      new("dPredModule", coef = numeric(ncol(from)),
          X = from, fac = chol(crossprod(from))))

setAs("dsparseModelMatrix", "predModule",
      function(from)
      new("sPredModule", coef = numeric(ncol(from)),
          X = from, fac = Cholesky(crossprod(from))))
