using Test
using Random
using StatsBase
using Distributions: Uniform
using OrdinaryDiffEqTsit5
using SciMLBase
using SymbolicIndexingInterface: SymbolCache
using Chemostats
using Chemostats: Offspring, ExactPartition, BinomialPartition, AllOrNothingPartition

# Integrator with constant state [V, P, R], to read parent values from
test_integrator(u = [4.0, 100.0, 1.0]) = init(ODEProblem((du, u, p, t) -> (du .= 0), u, (0., 1.)), Tsit5())
const iV, iP, iR = CartesianIndex(1), CartesianIndex(2), CartesianIndex(3)

daughters(offspring, int = test_integrator(); rng = Xoshiro(1)) = Chemostats.offspring_states(offspring, int, rng)

@testset "Number of offspring" begin
    @test_throws ArgumentError Offspring(iV => ExactPartition())     # no default number of offspring
    @test Offspring(3, iV => ExactPartition()).n == 3
    @test Offspring(iV => ExactPartition([0.2, 0.3, 0.5])).n == 3
    @test Offspring(iR => [1., 2., 3., 4.]).n == 4
    @test Offspring(0).n == 0

    @test_throws ArgumentError Offspring(-1)
    @test_throws ArgumentError Offspring(2, iV => [1., 2., 3.])
    @test_throws ArgumentError Offspring(iV => ExactPartition([0.5, 0.5]), iP => [1., 2., 3.])
    @test_throws ArgumentError Offspring(iV => 1.0, :P => 2.0)     # mixed indices and symbols
    @test_throws ArgumentError Offspring(1 => 1.0)                 # plain integers are ambiguous
end

@testset "Fractions" begin
    @test_throws ArgumentError ExactPartition([0.5, 0.6])
    @test_throws ArgumentError BinomialPartition([1.2, -0.2])
end

@testset "Rules" begin
    int = test_integrator()

    # Unlisted variables are inherited
    ds = daughters(Offspring(2, iV => ExactPartition()), int)
    @test [ d.u0 for d in ds ] == [ [2., 100., 1.], [2., 100., 1.] ]
    @test all(d -> isnothing(d.p), ds)

    ds = daughters(Offspring(iV => ExactPartition([0.25, 0.75])), int)
    @test [ d.u0[1] for d in ds ] == [1., 3.]

    ds = daughters(Offspring(iP => [5., 6., 7.], iR => 0.5), int)
    @test [ d.u0[2] for d in ds ] == [5., 6., 7.]
    @test all(d -> d.u0[3] == 0.5, ds)

    ds = daughters(Offspring(2, iR => Uniform(2, 3)), int)
    @test all(d -> 2 <= d.u0[3] <= 3, ds)
    @test ds[1].u0[3] != ds[2].u0[3]

    # Scalar states are set with `CartesianIndex()`
    scalar = init(ODEProblem((u, p, t) -> zero(u), 4.0, (0., 1.)), Tsit5())
    @test [ d.u0 for d in daughters(Offspring(2, CartesianIndex() => ExactPartition()), scalar) ] == [2., 2.]
    @test_throws ArgumentError daughters(Offspring(2, iV => 1.0), scalar)

    # Symbols give pairs, which `DECell` applies to the integrator
    ds = daughters(Offspring(:V => 1.0, :P => [3, 4]), int)
    @test ds[1].u0 == [ :V => 1.0, :P => 3 ]
    @test ds[2].u0 == [ :V => 1.0, :P => 4 ]

    # Rules for parameters set the parameters of the offspring
    sys = SymbolCache([:V], [:V_d], :t)
    pint = init(ODEProblem(ODEFunction((du, u, p, t) -> (du .= 0); sys), [4.0], (0., 1.), [2.0]), Tsit5())
    ds = daughters(Offspring(:V => ExactPartition(), :V_d => [1.0, 1.5]), pint)
    @test ds[1].u0 == [ :V => 2.0 ] && ds[1].p == [ :V_d => 1.0 ]
    @test ds[2].u0 == [ :V => 2.0 ] && ds[2].p == [ :V_d => 1.5 ]

    ds = daughters(Offspring(2, :V_d => ExactPartition()), pint)
    @test all(d -> isempty(d.u0) && d.p == [ :V_d => 1.0 ], ds)
end

@testset "Binomial partition" begin
    int = test_integrator()
    rng = Xoshiro(1)

    counts = [ [ d.u0[2] for d in daughters(Offspring(iP => BinomialPartition([0.3, 0.7])), int; rng) ] for _ in 1:2000 ]
    @test all(c -> sum(c) == 100 && all(>=(0), c) && all(isinteger, c), counts)
    @test mean(first.(counts)) ≈ 30 atol = 0.5

    counts = [ d.u0[2] for d in daughters(Offspring(4, iP => BinomialPartition()), int; rng) ]
    @test length(counts) == 4 && sum(counts) == 100
end

@testset "All-or-nothing partition" begin
    int = test_integrator()
    rng = Xoshiro(1)

    @test Offspring(iV => AllOrNothingPartition([0.2, 0.3, 0.5])).n == 3
    @test_throws ArgumentError AllOrNothingPartition([0.5, 0.6])

    volumes = [ [ d.u0[1] for d in daughters(Offspring(iV => AllOrNothingPartition([0.8, 0.2])), int; rng) ] for _ in 1:2000 ]
    @test all(v -> sort(v) == [0., 4.], volumes)
    @test mean(v -> v[1] == 4, volumes) ≈ 0.8 atol = 0.03
end
