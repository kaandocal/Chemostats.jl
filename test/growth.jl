using Test
using Random
using StatsBase
using Chemostats

include("models/exponential.jl")
include("models/sizecontrol.jl")
include("models/multitypemarkov.jl")

# Linear birth-death process: at each Exponential(1) event a cell either
# divides (probability p_divide) or dies without offspring. Growth rate is
# Λ = 2*p_divide - 1. Reuses f_exp/cb_exp from models/exponential.jl.
divide_bd(p_divide) = int -> rand() < p_divide ? map(_ -> (u0 = randexp(), p = (; p_divide)), 1:2) : nothing

function BirthDeathCell(p_divide; alg = Tsit5())
    prob = ODEProblem(f_exp, randexp(), (0., 0.), (; p_divide); callback = cb_exp)
    DECell(prob, alg, divide_bd(p_divide))
end

make_bd_population(K, p_divide; alg = Tsit5()) = [ BirthDeathCell(p_divide; alg) for _ in 1:K ]

rmse(Λs, Λ_gt) = sqrt(mean(abs2.(Λs .- Λ_gt)))

function rmse_Λ(make_chem, alg, Λ_gt, tmax, ensalg = Chemostats.default_ensalg(alg); niter = 10, kwargs...)
    Λs = map(1:niter) do _
        chem = make_chem()
        Chemostats.simulate!(chem, tmax, alg, ensalg; kwargs...)
        est_Λ(chem)
    end
    rmse(Λs, Λ_gt)
end

@testset "Yule process" begin
    Λ_gt = 1.0
    tmax = 100.0

    @testset "Strict(L=$L, clone=$clone)" for (L, tol) in [(20, 0.08), (100, 0.02), (250, 0.01)], clone in (true, false)
        Random.seed!(20260101 + L + 1000 * clone)
        make_chem = () -> Chemostat(make_population(L))
        @test rmse_Λ(make_chem, Chemostats.Strict(L, clone), Λ_gt, tmax) < tol
    end

    @testset "Lax(100) with $ensalg" for ensalg in (EnsembleSerial(), EnsembleThreads())
        Random.seed!(20260102)
        make_chem = () -> Chemostat([ ExponentialCell((; id = 1)) ])
        @test rmse_Λ(make_chem, Chemostats.Lax(250, 0.5), Λ_gt, tmax, ensalg; Nmax = 1e4) < 0.025
    end
end

@testset "Birth-death process (Strict clone on/off)" begin
    p_divide = 0.75
    Λ_gt = 2 * p_divide - 1
    tmax = 100.0

    @testset "Strict(L=$L, clone=$clone)" for (L, tol) in [(20, 0.15), (100, 0.05), (250, 0.03)], clone in (true, false)
        Random.seed!(20260301 + L + 1000 * clone)
        make_chem = () -> Chemostat(make_bd_population(L, p_divide))
        @test rmse_Λ(make_chem, Chemostats.Strict(L, clone), Λ_gt, tmax) < tol
    end
end

@testset "Size control (adder)" begin
    Λ_gt = 1.0
    tmax = 100.0

    @testset "Strict(L=$L, clone=$clone)" for (L, tol) in [(20, 0.08), (100, 0.02), (250, 0.01)], clone in (true, false)
        Random.seed!(20260201 + L + 1000 * clone)
        make_chem = () -> Chemostat(make_sizecontrol_population(L))
        @test rmse_Λ(make_chem, Chemostats.Strict(L, clone), Λ_gt, tmax) < tol
    end

    @testset "Lax(100) with $ensalg" for ensalg in (EnsembleSerial(), EnsembleThreads())
        Random.seed!(20260202)
        make_chem = () -> Chemostat([ SizeControlCell() ])
        @test rmse_Λ(make_chem, Chemostats.Lax(250, 0.5), Λ_gt, tmax, ensalg; Nmax = 1e4) < 0.025
    end
end

@testset "Multitype Markov" begin
    Random.seed!(20260101)
    prob = MultitypeMarkov([ 1., 2. ], [ 0.7 0.1; 0.3 0.9 ])
    Λ_gt = get_Λ_mtm(prob.p.model)
    alg = Tsit5()
    tmax = 100 / Λ_gt

    @testset "Strict(250, clone=$clone)" for clone in (true, false)
        Random.seed!(20260101)
        make_chem = () -> Chemostat([ DECell(prob, alg, divide_mtm) for _ in 1:250 ])
        @test rmse_Λ(make_chem, Chemostats.Strict(250, clone), Λ_gt, tmax) < 0.02 * Λ_gt
    end

    @testset "Lax(300) with $ensalg" for ensalg in (EnsembleSerial(), EnsembleThreads())
        Random.seed!(20260101)
        make_chem = () -> Chemostat([ DECell(prob, alg, divide_mtm) ])
        @test rmse_Λ(make_chem, Chemostats.Lax(300, 0.5 / Λ_gt), Λ_gt, tmax, ensalg; Nmax = 1e4) < 0.025 * Λ_gt
    end
end
