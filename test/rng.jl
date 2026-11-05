using Test
using Random
using SciMLBase
using Chemostats

include("models/exponential.jl")
include("models/jumpcell.jl")

# Exponential lifetimes; dies with probability `p_die`, otherwise divides into two.
# All randomness comes from the cell's RNG.
function divide_rng(int)
    rng = SciMLBase.get_rng(int)
    rand(rng) < int.p.p_die && return nothing
    map(_ -> (u0 = randexp(rng), p = int.p), 1:2)
end

function RNGCell(seed; p_die = 0.2)
    rng = Xoshiro(seed)
    prob = ODEProblem(f_exp, randexp(rng), (0., 0.), (; p_die); callback = cb_exp)
    DECell(prob, Tsit5(), divide_rng; rng)
end

@testset "Daughters and clones get independent streams of the same type" begin
    prob = ODEProblem(f_exp, 0.5, (0., 0.), (; p_die = 0.0); callback = cb_exp)
    cell = DECell(prob, Tsit5(), divide_rng; rng = MersenneTwister(1))
    Chemostats.init_cell!(cell)
    Chemostats.step!(cell, 0.1, nothing)
    clone = Chemostats.clone_cell(cell, Chemostats.get_curr_t(cell))
    Chemostats.step!(cell, 20.0, nothing)
    children = Chemostats.get_children(cell)

    @test length(children) == 2
    @test all(c -> typeof(c) == typeof(cell), [clone; children])
    @test allunique([ Chemostats.get_cell_rng(c) for c in [cell; clone; children] ])   # compares states
end

# Minimal cell for the queue's random selection
struct KeyCell
    t::Float64
    id::Int
    rng::Xoshiro
end

Chemostats.get_curr_t(c::KeyCell) = c.t
Chemostats.get_cell_rng(c::KeyCell) = c.rng
Chemostats.clone_cell(c::KeyCell, t) = KeyCell(t, c.id, Xoshiro(rand(c.rng, UInt64)))

# All at the same time, as at a tstop, so heap order alone can't break ties
key_cells(order) = [ KeyCell(1.0, i, Xoshiro(i)) for i in order ]

@testset "Random selection doesn't depend on heap layout" begin
    q1 = Chemostats.ThreadedQueue(key_cells(1:10))
    q2 = Chemostats.ThreadedQueue(key_cells(10:-1:1))
    @test [ v.id for v in Chemostats.iter_unsafe(q1) ] != [ v.id for v in Chemostats.iter_unsafe(q2) ]

    @test [ v.id for v in Chemostats._pop_random!(q1, 3) ] == [ v.id for v in Chemostats._pop_random!(q2, 3) ]
    @test first(Chemostats._clone_random!(q1, 2.0)).id == first(Chemostats._clone_random!(q2, 2.0)).id
end

# Population up to cell order, plus the full snapshot history
function fingerprint(chem)
    cells = sort([ (Chemostats.get_curr_t(c), copy(c.int.u)) for c in chem.pop ])
    snaps = [ (s.t, s.N, s.nsim, s.log_f) for s in chem.snaps ]
    cells, snaps
end

# The global RNG is scrambled first: results must depend only on `seed`
function run_population(ctor, K, alg, ensalg, tmax; seed = 1000)
    chem = Chemostat([ ctor(seed + i) for i in 1:K ])
    Random.seed!(rand(RandomDevice(), UInt64))
    Chemostats.simulate!(chem, tmax, alg, ensalg; saveat = 0.5:0.5:tmax)
    @test chem.status == SciMLBase.ReturnCode.Success
    fingerprint(chem)
end

const RNG_POP_CASES = [
    ("Thin", RNGCell, 20, Chemostats.Thin(0.3), 3.0),
    ("Forward", seed -> RNGCell(seed; p_die = 0.0), 10, Chemostats.Forward(10), 3.0),
    ("Lax", RNGCell, 40, Chemostats.Lax(40, 0.5), 5.0),
    ("Thin, jumps, remake", seed -> JumpCell(; k = 50.0, remake = true, rng = Xoshiro(seed)), 10, Chemostats.Thin(0.3), 2.5),
]

# Only meaningful with more than one thread, but must pass either way
@testset "simulate! is seed-reproducible across ensemble algorithms ($name)" for (name, ctor, K, alg, tmax) in RNG_POP_CASES
    ref = run_population(ctor, K, alg, EnsembleSerial(), tmax)
    @test isequal(run_population(ctor, K, alg, EnsembleThreads(), tmax), ref)
    @test !isequal(run_population(ctor, K, alg, EnsembleSerial(), tmax; seed = 2000), ref)
end

@testset "Strict is seed-reproducible" begin
    alg = Chemostats.Strict(20)
    ref = run_population(RNGCell, 20, alg, EnsembleSerial(), 4.0)
    @test isequal(run_population(RNGCell, 20, alg, EnsembleSerial(), 4.0), ref)
    @test !isequal(run_population(RNGCell, 20, alg, EnsembleSerial(), 4.0; seed = 2000), ref)
end
