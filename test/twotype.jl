using Test
using Chemostats 
using OrdinaryDiffEq

prob = Chemostats.Models.MultitypeMarkov([ 1., 2. ], [ 0.7 0.1; 0.3 0.9 ])
Λ_gt = Chemostats.Models.get_Λ(prob)

###

# niter = 1000
# @testset "Extinction probability (thin, δ = $δ)" for δ in [ 0., 0.5, 0.9, 0.99 ]
#     tmax = 4 / (1 - δ)
#     NN = map(1:niter) do i
#         chem = Chemostats.Chemostat([ Chemostats.Models.sample_cell(prob) ])
#         Chemostats.simulate!(chem, tmax, Chemostats.Thin(δ * Λ_gt), ensalg) 
#         chem.saved[end].N
#     end

#     @test count(iszero, NN) / niter ≈ δ atol=0.01
# end 

### 

tmax = 50.
tt = 0:0.1:tmax

@testset "System size (strict, L=$L)" for L in [ 1, 10 ]
    chem = Chemostats.Chemostat([ Chemostats.Models.sample_cell(prob) ])
    if ensalg isa EnsembleThreads 
        @test_throws "does not support parallelisation" Chemostats.simulate!(chem, tmax, Chemostats.Strict(L), ensalg; saveat=tt)
    else 
        Chemostats.simulate!(chem, tmax, Chemostats.Strict(L), ensalg; saveat=tt)

        Ns = [ snap.N for snap in chem.saved[2:end] ]
        @test all(Ns .== L)
    end
end

@testset "System size (forward, L=$L)" for L in [ 1, 10 ]
    chem = Chemostats.Chemostat([ Chemostats.Models.sample_cell(prob) for i in 1:L ])
    Chemostats.simulate!(chem, tmax, Chemostats.Forward(L), ensalg; saveat=tt)

    Ns = [ snap.N for snap in chem.saved ]
    @test all(Ns .== L)
end

### 

@testset "Λ estimator (strict, L=10k)" begin
    tmax = 400 / Λ_gt 
    niter = 10

    ΛΛ = map(1:niter) do i
        chem = Chemostats.Chemostat([ Chemostats.Models.sample_cell(prob) for i in 1:10000 ])
        Chemostats.simulate!(chem, tmax, Chemostats.Strict(10000), ensalg) 
        Chemostats.est_Λ(chem)
    end

    @test count(abs.(ΛΛ .- Λ_gt) .< 0.01 * Λ_gt) / niter >= 0.7
end 

@testset "Λ estimator (lax, L=10k)" begin
    tmax = 400 / Λ_gt 
    niter = 10

    ΛΛ = map(1:niter) do i
        chem = Chemostats.Chemostat([ Chemostats.Models.sample_cell(prob) for i in 1:10000 ])
        Chemostats.simulate!(chem, tmax, Chemostats.Lax(10000, 2), ensalg) 
        Chemostats.est_Λ(chem)
    end

    @test count(abs.(ΛΛ .- Λ_gt) .< 0.01 * Λ_gt) / niter >= 0.7
end 