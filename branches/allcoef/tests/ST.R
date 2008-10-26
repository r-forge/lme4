## Elementary tests of the ST class

library(lme4)
tst <- function(ST, Zt, pars)
{
    stopifnot(is(ST, "ST"), is(Zt, "dgCMatrix"))
    .Call("ST_initialize", ST, Zt)
    print(str(ST))
    print(.Call("ST_getPars", ST))
    print(.Call("ST_bounds", ST))
    .Call("ST_setPars", ST, pars)
    print(.Call("ST_Tmatrix", ST))
    print(Lambda <- .Call("ST_Lambda", ST))
    print(t(Lambda) %*% Zt)
    print(str(A <- .Call("ST_create_A", ST, Zt)))
    .Call("ST_update_A", ST, Zt, A)
    print(str(A))
}
tst(new("ST", ST = list(matrix(1)), Gp = c(0L, 6L)),
    as(gl(6,5), "sparseMatrix"), 0.8)
tst(new("ST", ST = list(matrix(1), matrix(1)), Gp = c(0L, 3L, 6L)),
    rbind2(as(gl(3, 9), "sparseMatrix"), as(gl(3, 3, 27), "sparseMatrix")),
    c(0.8,0.6))
Zt1 <- Zt2 <- as(gl(3, 5), "sparseMatrix")
Zt2@x <- as.double(rep(1:5, 3))
st3 <- new("ST", ST = list(diag(nrow = 2)), Gp = c(0L, 6L))
tst(st3, rbind2(Zt1, Zt2), c(0.8, 0.2, -0.1))
