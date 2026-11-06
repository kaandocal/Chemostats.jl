using Test
using Random
using Statistics
using Distributions: Uniform
using OrdinaryDiffEqTsit5
using SciMLBase
using Chemostats
using Chemostats: Offspring, ExactPartition, BinomialPartition, AllOrNothingPartition

# Integrator with constant state [V, P, R], to read parent values from
test_integrator(u = [4.0, 100.0, 1.0]) = init(ODEProblem((du, u, p, t) -> (du .= 0), u, (0., 1.)), Tsit5())

daughters(offspring, int = test_integrator(); rng = Xoshiro(1)) = Chemostats.offspring_states(offspring, int, rng)

@testset "Number of offspring" begin
    @test Offspring(1 => ExactPartition()).n == 2
    @test Offspring(3, 1 => ExactPartition()).n == 3
    @test Offspring(1 => ExactPartition([0.2, 0.3, 0.5])).n == 3
    @test Offspring(3 => [1., 2., 3., 4.]).n == 4
    @test Offspring(0).n == 0

    @test_throws ArgumentError Offspring(-1)
    @test_throws ArgumentError Offspring(2, 1 => [1., 2., 3.])
    @test_throws ArgumentError Offspring(1 => ExactPartition([0.5, 0.5]), 2 => [1., 2., 3.])
    @test_throws ArgumentError Offspring(1 => 1.0, :P => 2.0)     # mixed indices and symbols
end

@testset "Fractions" begin
    @test_throws ArgumentError ExactPartition([0.5, 0.6])
    @test_throws ArgumentError BinomialPartition([1.2, -0.2])
end

@testset "Rules" begin
    int = test_integrator()

    # Unlisted variables are inherited
    ds = daughters(Offspring(1 => ExactPartition()), int)
    @test [ d.u0 for d in ds ] == [ [2., 100., 1.], [2., 100., 1.] ]
    @test all(d -> isnothing(d.p), ds)

    ds = daughters(Offspring(1 => ExactPartition([0.25, 0.75])), int)
    @test [ d.u0[1] for d in ds ] == [1., 3.]

    ds = daughters(Offspring(2 => [5., 6., 7.], 3 => 0.5), int)
    @test [ d.u0[2] for d in ds ] == [5., 6., 7.]
    @test all(d -> d.u0[3] == 0.5, ds)

    ds = daughters(Offspring(3 => Uniform(2, 3)), int)
    @test all(d -> 2 <= d.u0[3] <= 3, ds)
    @test ds[1].u0[3] != ds[2].u0[3]

    # Symbols give pairs, which `DECell` applies to the integrator
    ds = daughters(Offspring(:V => 1.0, :P => [3, 4]), int)
    @test ds[1].u0 == [ :V => 1.0, :P => 3 ]
    @test ds[2].u0 == [ :V => 1.0, :P => 4 ]
end

@testset "Binomial partition" begin
    int = test_integrator()
    rng = Xoshiro(1)

    counts = [ [ d.u0[2] for d in daughters(Offspring(2 => BinomialPartition([0.3, 0.7])), int; rng) ] for _ in 1:2000 ]
    @test all(c -> sum(c) == 100 && all(>=(0), c) && all(isinteger, c), counts)
    @test mean(first.(counts)) ≈ 30 atol = 0.5

    counts = [ d.u0[2] for d in daughters(Offspring(4, 2 => BinomialPartition()), int; rng) ]
    @test length(counts) == 4 && sum(counts) == 100
end

@testset "All-or-nothing partition" begin
    int = test_integrator()
    rng = Xoshiro(1)

    @test Offspring(1 => AllOrNothingPartition([0.2, 0.3, 0.5])).n == 3
    @test_throws ArgumentError AllOrNothingPartition([0.5, 0.6])

    volumes = [ [ d.u0[1] for d in daughters(Offspring(1 => AllOrNothingPartition([0.8, 0.2])), int; rng) ] for _ in 1:2000 ]
    @test all(v -> sort(v) == [0., 4.], volumes)
    @test mean(v -> v[1] == 4, volumes) ≈ 0.8 atol = 0.03
end

@testset "Division of a DECell" begin
    # The volume grows exponentially, and the cell divides once it has doubled
    grow!(du, u, p, t) = (du[1] = u[1])
    cb = Chemostats.DivideCallback((u, t, int) -> u[1] - 2.0)
    prob = ODEProblem(grow!, [1.0], (0., Inf); callback = cb)

    function children(divide; seed = 1)
        cell = DECell(prob, Tsit5(), divide; rng = Xoshiro(seed))
        Chemostats.init_cell!(cell)
        Chemostats.step!(cell, 10.0, nothing)
        Chemostats.get_children(cell)
    end
    volume(cell) = cell.int.u[1]

    kids = children(Offspring(1 => ExactPartition()))
    @test length(kids) == 2
    @test all(k -> isapprox(volume(k), 1; atol = 1e-6), kids)

    @test length(children(int -> Offspring(3, 1 => ExactPartition()))) == 3
    @test isempty(children(int -> nothing))
    @test isempty(children(Offspring(0)))

    # Random rules use the RNG of the dividing cell
    sample(seed) = volume.(children(Offspring(1 => Uniform(0.5, 1.5)); seed))
    @test sample(1) == sample(1)
    @test sample(1) != sample(2)
end
