using Test
using StatsBase
using Chemostats

model = model_vol_50
env = nothing 

Λ_gt = get_Λ(model)

### 

@testset "Λ estimator (thin, δ = 0.9)" begin
    tmax = 50.
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
    tmax = 600.
    niter = 10

    ΛΛ = map(1:niter) do i
        chem = Chemostats.Chemostat([ draw_cell(model, env) for i in 1:10000 ], model, env)
        Chemostats.simulate!(chem, tmax, Chemostats.Thin(Λ_gt)) 
        Chemostats.est_Λ(chem)
    end

    @test mean(ΛΛ) ≈ Λ_gt atol=0.001
end