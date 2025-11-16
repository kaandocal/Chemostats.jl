using Test
using Distributions
using Chemostats 
using OrdinaryDiffEq

f(u, p, t) = u * p.λ
get_Λ(p) = p.λ

p = (; Vd = 2., λ = log(2))
u0 = p.Vd / 2

divide_condition(u, t, int) = u - int.p.Vd
cb = Chemostats.DivideCallback(divide_condition)

prob = ODEProblem(f, u0, (0., 0.), p; callback=cb)
Λ_gt = get_Λ(p)

function Chemostats.get_offspring(int, cell::Chemostats.DECell; save_lineages=false)
    anc = if save_lineages
        cell 
    else 
        missing
    end 

    map(1:2) do _
        u0 = cell.int.u / 2
        prob = remake(cell.int.sol.prob; u0)
        Chemostats.DECell(anc, prob; t0=Chemostats.get_curr_t(cell))
    end
end 

model = prob 
env = nothing 

if !isdefined(Main, :ensalg)
    ensalg = Chemostats.EnsembleSerial()
end

###

niter = 1000
@testset "Extinction probability (thin, δ = $δ)" for δ in [ 0.5, 0.9, 0.99 ]
    tmax = 200. * δ^3
    NN = map(1:niter) do i
        chem = Chemostats.Chemostat([ Chemostats.DECell(prob) ], model, env)
        Chemostats.simulate!(chem, tmax, Chemostats.Thin(δ * Λ_gt), ensalg) 
        chem.saved[end].N
    end

    @test count(iszero, NN) / niter ≈ exp(δ * Λ_gt) - 1 atol=0.02
end 

### 

@testset "System size (direct)" begin 
    tmax = 11.5

    chem = Chemostats.Chemostat([ Chemostats.DECell(prob) ], model, env)
    Chemostats.simulate!(chem, tmax, Chemostats.Direct(), ensalg) 
    N = chem.saved[end].N

    @test N == 2 ^ floor(Int, tmax)
end 

### 

tmax = 50.
tt = 0:0.1:tmax

@testset "System size (strict, L=$L)" for L in [ 1, 10 ]
    chem = Chemostats.Chemostat([ Chemostats.DECell(prob) ], model, env)
    Chemostats.simulate!(chem, tmax, Chemostats.Strict(L), ensalg; saveat=tt)

    Ns = [ snap.N for snap in chem.saved[2:end] ]
    @test all(Ns .== L)
end

@testset "System size (forward, L=$L)" for L in [ 1, 10 ]
    chem = Chemostats.Chemostat([ Chemostats.DECell(prob) for i in 1:L ], model, env)
    Chemostats.simulate!(chem, tmax, Chemostats.Forward(L), ensalg; saveat=tt)

    Ns = [ snap.N for snap in chem.saved ]
    @test all(Ns .== L)
end

### 

@testset "Λ estimator (thin, δ = 0.9)" begin
    tmax = 100
    niter = 100

    ΛΛ = map(1:niter) do i
        chem = Chemostats.Chemostat([ Chemostats.DECell(prob) for i in 1:100 ], model, env)
        Chemostats.simulate!(chem, tmax, Chemostats.Thin(0.9 * Λ_gt), ensalg) 
        Chemostats.est_Λ(chem)
    end

    @test count(abs.(ΛΛ .- Λ_gt) .< 0.01) / niter >= 0.7
end 

###

@testset "Λ estimator (thin, δ = 0.99)" begin
    tmax = 200.
    niter = 100

    ΛΛ = map(1:niter) do i
        chem = Chemostats.Chemostat([ Chemostats.DECell(prob) for i in 1:1000 ], model, env)
        Chemostats.simulate!(chem, tmax, Chemostats.Thin(0.99 * Λ_gt), ensalg) 
        Chemostats.est_Λ(chem)
    end

    @test count(abs.(ΛΛ .- Λ_gt) .< 0.005) / niter >= 0.75
end

### 

@testset "Λ estimator (thin, δ = 1)" begin 
    tmax = 1000.
    niter = 10

    ΛΛ = map(1:niter) do i
        chem = Chemostats.Chemostat([ Chemostats.DECell(prob) for i in 1:10000 ], model, env)
        Chemostats.simulate!(chem, tmax, Chemostats.Thin(Λ_gt), ensalg) 
        Chemostats.est_Λ(chem)
    end

    @test mean(ΛΛ) ≈ Λ_gt atol=0.001
end

###

tmax = 100.5

@testset "Reaction counts (forward, L=$L)" for L in [ 1, 10, 100 ]
    chem = Chemostats.Chemostat([ Chemostats.DECell(prob) for i in 1:L ], model, env)
    Chemostats.simulate!(chem, tmax, Chemostats.Forward(L), ensalg) 
    
    n = chem.saved[end].log_f / log(1+1/L)

    @test n ≈ floor(Int, tmax) * L
end

