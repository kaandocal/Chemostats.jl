using Test
using Random
using Unzip
using Distributions
using SpecialFunctions

using Revise
include("../src/Chemostats.jl")
Cell = Chemostats.Cell
kill_cell! = Chemostats.kill_cell!
function Chemostats.simulate_cell!(args...; kwargs...) 
    simulate_cell!(args...; kwargs...)
end 

divide! = Chemostats.divide!
CellState = Chemostats.CellState 

function Chemostats.get_offspring(args...; kwargs...) 
    get_offspring(args...; kwargs...)
end 

###

test_poisson(samples, μ_pred; kwargs...) = test_dist(samples, Poisson(μ_pred); kwargs...)
function test_dist(samples, dist_pred; α=0.05)
    zm, zv = zscore_meanvar(samples, dist_pred)

    tol = -quantile(Normal(), α/2)

    @test abs(zm) <= tol
    @test abs(zv) <= tol
end 

function zscore_meanvar(samples, dist)
    n = length(samples)
    m_gt = mean(dist)
    v_gt = var(dist)

    m_emp = mean(samples)
    v_emp = mean(abs2.(samples .- m_gt))

    zm = (m_emp - m_gt) * sqrt(n) / std(dist)
    zv = (v_emp - v_gt) * sqrt(n) / sqrt((3 + kurtosis(dist)) * var(dist)^2)

    zm, zv
end 

###

include("../models/matrix.jl")

model = model_yule
env = nothing 

Λ_gt = get_Λ(model)

###

niter = 1000
@testset "Extinction probability (thin, δ = $δ)" for δ in [ 0., 0.5, 0.9, 0.99 ]
    tmax = 4 / (1 - δ)
    NN = map(1:niter) do i
        chem = Chemostats.Chemostat([ draw_cell(model, env) ], model, env)
        Chemostats.simulate!(chem, tmax, Chemostats.Thin(δ * Λ_gt)) 
        chem.saved[end].N
    end

    @test count(iszero, NN) / niter ≈ δ atol=0.01
end 

### 

@testset "System size (thin, δ = $δ)" for δ in [ 0., 0.5, 0.9 ]#, 0.99 ]
    tmax = 5 / (1 - δ)
    niter = 1000 / (1 - δ)
    NN = map(1:niter) do i
        chem = Chemostats.Chemostat([ draw_cell(model, env) ], model, env)
        Chemostats.simulate!(chem, tmax, Chemostats.Thin(δ * Λ_gt)) 
        chem.saved[end].N
    end

    filter!(!iszero, NN)

    mean_pred = exp((1 - δ) * Λ_gt * tmax) / (1 - δ)

    dist_pred = Geometric(1 / (mean_pred + 1))

    test_dist(NN, dist_pred)
end 

### 

tmax = 50.
tt = 0:0.1:tmax

@testset "System size (strict, L=$L)" for L in [ 1, 10 ]
    chem = Chemostats.Chemostat([ draw_cell(model, env) ], model, env)
    Chemostats.simulate!(chem, tmax, Chemostats.Strict(L); saveat=tt)

    Ns = [ snap.N for snap in chem.saved[2:end] ]
    @test all(Ns .== L)
end

@testset "System size (forward, L=$L)" for L in [ 1, 10 ]
    chem = Chemostats.Chemostat([ draw_cell(model, env) for i in 1:L ], model, env)
    Chemostats.simulate!(chem, tmax, Chemostats.Forward(L); saveat=tt)

    Ns = [ snap.N for snap in chem.saved ]
    @test all(Ns .== L)
end

### 

@testset "Λ estimator (thin, δ = 0.9)" begin
    tmax = 50
    niter = 100

    ΛΛ = map(1:niter) do i
        chem = Chemostats.Chemostat([ draw_cell(model, env) for i in 1:100 ], model, env)
        Chemostats.simulate!(chem, tmax, Chemostats.Thin(0.9 * Λ_gt)) 
        Chemostats.est_Λ(chem)
    end

    @test count(abs.(ΛΛ .- Λ_gt) .< 0.01) / niter >= 0.7
end 

###

@testset "Λ estimator (thin, δ = 0.99)" begin
    tmax = 100.
    niter = 100

    ΛΛ = map(1:niter) do i
        chem = Chemostats.Chemostat([ draw_cell(model, env) for i in 1:1000 ], model, env)
        Chemostats.simulate!(chem, tmax, Chemostats.Thin(0.99 * Λ_gt)) 
        Chemostats.est_Λ(chem)
    end

    @test count(abs.(ΛΛ .- Λ_gt) .< 0.005) / niter >= 0.75
end

### 

@testset "Λ estimator (thin, δ = 1)" begin 
    tmax = 30.
    niter = 10

    ΛΛ = map(1:niter) do i
        chem = Chemostats.Chemostat([ draw_cell(model, env) for i in 1:50000 ], model, env)
        Chemostats.simulate!(chem, tmax, Chemostats.Thin(Λ_gt)) 
        Chemostats.est_Λ(chem)
    end

    @test mean(ΛΛ) ≈ Λ_gt atol=0.001
end

### 

tmax = 100.
niter = 1000

@testset "Reaction counts (strict, L=$L)" for L in [ 1, 10, 100 ]
    nn = map(1:niter) do i
        chem = Chemostats.Chemostat([ draw_cell(model, env) for i in 1:L ], model, env)
        Chemostats.simulate!(chem, tmax, Chemostats.Strict(L)) 
        chem.saved[end].log_f / log(1+1/L)
    end

    test_poisson(nn, Λ_gt * tmax * L)
end

@testset "Reaction counts (forward, L=$L)" for L in [ 1, 10, 100 ]
    nn = map(1:niter) do i
        chem = Chemostats.Chemostat([ draw_cell(model, env) for i in 1:L ], model, env)
        Chemostats.simulate!(chem, tmax, Chemostats.Forward(L)) 
        chem.saved[end].log_f / log(1+1/L)
    end

    test_poisson(nn, Λ_gt * tmax * L)
end

