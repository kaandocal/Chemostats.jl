using Test
using Chemostats

include("models/exponential.jl")

ENSALGS = (EnsembleSerial(), EnsembleThreads())

# Test the scheduler behaves under edge conditions (explosion, extinction, errors)
@testset "Nmax abort" for ensalg in ENSALGS
    K = 50
    chem = Chemostat(make_population(K))
    alg = Chemostats.Direct()

    int = Chemostats.PopIntegrator(chem, alg, ensalg)
    Chemostats.init!(alg, int)
    Chemostats.simulate!(int, 20.0, ensalg; Nmax = K)

    @test int.retcode == Chemostats.ReturnCode.MaxIters
    @test length(chem.pop) > 0
end

@testset "Population extinction" begin
    EXTINCTION_ALGS = (
        Chemostats.Direct(),
        Chemostats.Thin(0.5),
        Chemostats.Forward(15),
        Chemostats.Lax(15, 0.5),
    )

    for alg in EXTINCTION_ALGS, ensalg in ENSALGS
        @testset "$(typeof(alg).name.name) / $ensalg" begin
            K = 15
            chem = Chemostat(make_dying_population(K))
            tmax = 20.0

            Chemostats.simulate!(chem, tmax, alg, ensalg; Nmax = 10_000)
            snap = chem.snaps[end]

            @test snap.t == tmax
            @test snap.N == 0
            @test length(chem.pop) == 0
        end
    end
end

@testset "throw_on_error handling" begin
    K = 10
    tmax = 15.0

    @testset "throw_on_error = true ($ensalg)" for ensalg in ENSALGS
        chem = Chemostat(make_throwing_population(K))
        caught = nothing
        try
            Chemostats.simulate!(chem, tmax, Chemostats.Direct(), ensalg; throw_on_error = true)
        catch e
            caught = e
        end

        @test !isnothing(caught)
        @test has_cause(caught, CellError)
        @test has_cause(caught, CellException)   # the divide error is wrapped, not raw
    end

    @testset "throw_on_error = false ($ensalg)" for ensalg in ENSALGS
        chem = Chemostat(make_throwing_population(K))

        Chemostats.simulate!(chem, tmax, Chemostats.Direct(), ensalg)

        @test length(chem.pop) == 0

        @test chem.snaps[end].t == tmax
    end

    @testset "throw_on_error = false retcode ($ensalg)" for ensalg in ENSALGS
        chem = Chemostat(make_throwing_population(K))
        alg = Chemostats.Direct()
        int = Chemostats.PopIntegrator(chem, alg, ensalg)
        Chemostats.init!(alg, int)

        Chemostats.simulate!(int, tmax, ensalg; throw_on_error = false)

        @test int.retcode == Chemostats.ReturnCode.Success
        @test length(chem.pop) == 0
    end

    @testset "throw_on_error=false throws on non-CellException ($ensalg)" for ensalg in ENSALGS
        chem = Chemostat([ BuggyCell(Chemostats.CellState.Newborn, 0.0) ])
        caught = nothing
        try
            Chemostats.simulate!(chem, 1.0, Chemostats.Direct(), ensalg; throw_on_error = false)
        catch e
            caught = e
        end

        @test !isnothing(caught)
        @test has_cause(caught, ErrorException)
        @test !has_cause(caught, CellException)
    end
end

@testset "Mixed offspring counts" begin
    K = 20
    tmax = 3.0

    MIXED_ALGS = (
        Chemostats.Direct(),
        Chemostats.Thin(0.1),
        Chemostats.Strict(K),
        Chemostats.Lax(K, 0.5),
    )

    for alg in MIXED_ALGS, ensalg in ENSALGS
        alg isa Chemostats.Strict && ensalg isa EnsembleThreads && continue

        @testset "$(typeof(alg).name.name) / $ensalg" begin
            chem = Chemostat(make_stochastic_population(K))
            Chemostats.simulate!(chem, tmax, alg, ensalg; Nmax = 10_000)
            snap = chem.snaps[end]

            @test snap.t == tmax
            @test snap.N > 0
            @test length(chem.pop) == snap.N
            @test snap.nsim >= K

            # Neither division nor cloning ever introduces a new id
            @test all(cell -> cell_params(cell).id in 1:K, chem.pop)

            if alg isa Chemostats.Strict
                @test snap.N == K
            end
        end
    end
end
