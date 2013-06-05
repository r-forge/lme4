using DataFrames
using Base.LinAlg.BLAS: trmm, gemv!

type MixedModel{Ti<:Union(Int32,Int64)}
    X::ModelMatrix                      # (dense) model matrix
    Xs::Vector{Matrix{Float64}}         # X_1,X_2,...,X_k
    beta::Vector{Float64}
    inds::Matrix{Ti}                    # n by k
    lambda::Vector{Matrix{Float64}}     # k lower triangular mats
    u::Vector{Matrix{Float64}}
    y::Vector{Float64}
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
