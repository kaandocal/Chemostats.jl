using Test
using Chemostats 

include("models/multitypemarkov.jl")

prob = MultitypeMarkov([ 1., 2. ], [ 0.7 0.1; 0.3 0.9 ])
Λ_gt = get_Λ_mtm(prob.p.model)

@testset "Strict(100)" begin
    tmax = 100 / Λ_gt 
    niter = 10

    ΛΛ = map(1:niter) do i
        chem = Chemostat([ DECell(prob, divide_mtm) for i in 1:100 ])
        Chemostats.simulate!(chem, tmax, Chemostats.Strict(100)) 
        est_Λ(chem)
    end

    @test sqrt(mean(abs2.(ΛΛ .- Λ_gt))) < 0.02 * Λ_gt
end

for ensalg in [ EnsembleSerial(), EnsembleThreads() ]
    @testset "Lax(100) with $ensalg" begin
        tmax = 100 / Λ_gt 
        niter = 10

        ΛΛ = map(1:niter) do i
            chem = Chemostat([ DECell(prob, divide_mtm) ])
            Chemostats.simulate!(chem, tmax, Chemostats.Lax(100, 0.5 / Λ_gt), ensalg; Nmax=1e4) 
            est_Λ(chem)
        end

        @test sqrt(mean(abs2.(ΛΛ .- Λ_gt))) < 0.04 * Λ_gt
    end 
end 
