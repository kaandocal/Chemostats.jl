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

include("../models/vol.jl")

model = model_det
env = nothing 

Λ_gt = get_Λ(model)

###

niter = 1000
@testset "Extinction probability (thin, δ = $δ)" for δ in [ 0.5, 0.9, 0.99 ]#[ 0., 0.5, 0.9 ]#, 0.99 ]
    tmax = 120. * δ^2 
    NN = map(1:niter) do i
        chem = Chemostats.Chemostat([ draw_cell(model, env) ], model, env)
        Chemostats.simulate!(chem, tmax, Chemostats.Thin(δ * Λ_gt)) 
        chem.saved[end].N
    end

    @test count(iszero, NN) / niter ≈ exp(δ * Λ_gt) - 1 atol=0.02
end 

### 

@testset "System size (thin, δ = 0)" begin 
    tmax = 11.5

    chem = Chemostats.Chemostat([ draw_cell(model, env) ], model, env)
    Chemostats.simulate!(chem, tmax, Chemostats.Thin(0)) 
    N = chem.saved[end].N

    @test N == 2 ^ floor(Int, tmax)
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
    Chemostats.simulate!(chem, tmax, Chemostats.Forward(); saveat=tt)

    Ns = [ snap.N for snap in chem.saved ]
    @test all(Ns .== L)
end

### 

@testset "Λ estimator (thin, δ = 0.9)" begin
    tmax = 100
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
    tmax = 200.
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
    tmax = 1000.
    niter = 10

    ΛΛ = map(1:niter) do i
        chem = Chemostats.Chemostat([ draw_cell(model, env) for i in 1:10000 ], model, env)
        Chemostats.simulate!(chem, tmax, Chemostats.Thin(Λ_gt)) 
        Chemostats.est_Λ(chem)
    end

    @test mean(ΛΛ) ≈ Λ_gt atol=0.001
end

###

tmax = 100.5

@testset "Reaction counts (forward, L=$L)" for L in [ 1, 10, 100 ]
    chem = Chemostats.Chemostat([ draw_cell(model, env) for i in 1:L ], model, env)
    Chemostats.simulate!(chem, tmax, Chemostats.Forward()) 
    
    n = chem.saved[end].log_f / log(1+1/L)

    @test n ≈ floor(Int, tmax) * L
end

