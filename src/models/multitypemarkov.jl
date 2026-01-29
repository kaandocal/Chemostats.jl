using Random
using StatsBase
using OrdinaryDiffEq
using LinearAlgebra
using ArgCheck
using UnPack
using Chemostats

f_mtm(u, p, t) = -1.
cb_mtm = ContinuousCallback((u, t, int) -> u, terminate!, interp_points=0)

function divide_mtm(int)
    @unpack s, model = int.p
    @unpack λ, T = model 

    if rand() > sum(T[:,s])
        return nothing 
    end 
    
    w = ProbabilityWeights(T[:,s])

    map(1:2) do _
        s_new = sample(w)
        u0_new = randexp() / λ[s_new]
        (; u0=u0_new, p=(; s=s_new, model))
    end
end 

function MultitypeMarkov(λ::AbstractVector, T::AbstractMatrix)
    @argcheck size(T) == (length(λ), length(λ))
    @argcheck all(sum(T; dims=1) .== 1)

    s = sample(1:length(λ))
    u0 = randexp() / λ[s]
    p = (; s, model = (; λ, T))

    ODEProblem(f_mtm, u0, (0., 0.), p; callback=cb_mtm)
end 

Yule(λ = 1.) = MultitypeMarkov([ λ ], [ 1;; ])

function sample_cell_mtm(prob)
    @unpack model = prob.p
    @unpack λ, T = model 

    s = sample(1:length(λ))
    u0 = randexp() / λ[s]

    prob_ = remake(prob; u0=u0, p=(; s, model))
    Chemostats.DECell(prob_, divide_mtm)
end 

function get_Λ_mtm(model::NamedTuple) 
    Q = 2 .* model.T .* model.λ - diagm(model.λ)
    last(eigvals(Q))
end

