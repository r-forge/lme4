using DataFrames
using Base.LinAlg.BLAS: trmm, gemm!
using Base.LinAlg.CHOLMOD: CholmodFactor, chm_factorize_p!

type MixedModel{Ti<:Union(Int32,Int64)}
    L::CholmodFactor{Float64,Ti}
    LambdatZt::CholmodSparse{Float64,Ti}
    X::ModelMatrix                      # fixed-effects model matrix
    Xs::Vector{Matrix{Float64}}         # X_1,X_2,...,X_k
    beta::Vector{Float64}
    inds::Matrix{Ti}                    # n by k
    lambda::Vector{Matrix{Float64}}     # k lower triangular mats
    u::Vector{Matrix{Float64}}
    y::Vector{Float64}
end

function updateL(m::MixedModel, theta::Vector{Float64})
    n,k = size(m.inds)
    LamZ = m.LambdatZt
    N = reshape(LamZ.nzvals, (div(length(LamZ),n),n))
    lambda = m.lambda; Xs = m.Xs
    tpos = 1; roff = 0                  # position in theta, row offset
    for i in 1:k
        T = lambda[i]
        p = size(T,1)                   # size of i'th template matrix
        for j in 1:p, i in j:p          # fill lower triangle from theta
            T[i,j] = theta[tpos]; tpos += 1
            if i == j && T[i,j] < 0. error("Negative diagonal element in T") end
        end
        gemm!('T','T',1.0,T,Xs[i],0.0,sub(nzmat,roff+(1:p),1:n))
        roff += p
    end
    chm_factorize_p!(m.L,m.LambdaZt,1.)
    m
end
        
function linpred(m::MixedModel)
    lp = m.X.m * m.beta                 # initialize value to X*beta
    lm = m.lambda
    Xs = m.Xs
    for i in 1:length(Xs)               # iterate over r.e. terms
        X = Xs[i]
        gemv!('N', 1.,
              trmm('R','L','N','N',1.,lm[i],X) .* m.u[i][m.inds[:,i],:],
              ones(size(X,2)), 1., lp)
    end
    lp
end

