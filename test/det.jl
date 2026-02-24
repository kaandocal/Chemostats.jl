using Test
using Distributions
using Chemostats 

include("models/deterministic.jl")

Λ_gt = 1.4
prob = Deterministic(Λ_gt)
alg = Tsit5()

for ensalg in [ EnsembleSerial(), EnsembleThreads() ]
    @testset "Direct() with $ensalg" begin 
        tmax = 5. / Λ_gt

        chem = Chemostat([ DECell(prob, alg, divide_det) ])
        Chemostats.simulate!(chem, tmax, Chemostats.Direct(), ensalg) 
        snap = chem.snaps[end]

        @test snap.t == tmax
        @test snap.N == 2 ^ floor(Int, tmax * Λ_gt)
        @test snap.nsim == 2 * snap.N - 1
        @test Chemostats.est_N(snap) == snap.N
    end 
end 

### 

tmax = 20.5 

@testset "Strict($L)" for L in [ 1, 10 ]
    chem = Chemostat([ DECell(prob, alg, divide_det) for i in 1:L ])

    @test_throws "Strict does not support parallelisation" Chemostats.simulate!(chem, tmax, Chemostats.Strict(L), EnsembleThreads())

    Chemostats.simulate!(chem, tmax, Chemostats.Strict(L))
    snap = chem.snaps[end]

    @test snap.t == tmax
    @test snap.N == L

    if L == 1
        Chemostats.est_logN(snap) ≈ floor(Int, tmax * Λ_gt) * log(1 + 1/L) * L + log(L)
    end 
end

for ensalg in [ EnsembleSerial(), EnsembleThreads() ]
    @testset "Forward($L) with $ensalg" for L in [ 1, 10 ]
        chem = Chemostat([ DECell(prob, alg, divide_det) for i in 1:L ])
        
        @test_throws DimensionMismatch Chemostats.simulate!(chem, tmax, Chemostats.Forward(L+1), ensalg)

        Chemostats.simulate!(chem, tmax, Chemostats.Forward(L))
        snap = chem.snaps[end]

        @test snap.t == tmax
        @test snap.N == L
        @test isnan(snap.log_f)
        @test snap.nsim == ceil(Int, tmax * Λ_gt) * L
    end
end 
