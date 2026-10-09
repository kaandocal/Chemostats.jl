using Test
using Random
using SciMLBase
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using Chemostats

include("models/jumpcell.jl")
include("models/flipcell.jl")

function run_to_division(cell)
    Chemostats.init_cell!(cell)
    Chemostats.step!(cell, 10.0, nothing)
    cell
end

run_to_division(prob::SciMLBase.AbstractDEProblem) = run_to_division(DECell(prob, Tsit5(), divide_jumpcell))

hybrid_prob(callback) = HybridProblem(jump_rn, [:P => 0, :V => 1.0], (0., Inf), [:k => 1.0, :V_d => 2.0]; callback)

@testset "DivideCallback with MTK equations" begin
    ref = run_to_division(JumpCell())

    for cond in (jump_rn.V ~ jump_rn.V_d, [ jump_rn.V ~ jump_rn.V_d ])
        cb_eq = Chemostats.DivideCallback(cond; interp_points=0, save_positions=(true, false))
        cell = run_to_division(hybrid_prob(cb_eq))

        @test Chemostats.get_state(cell) == Chemostats.CellState.EndOfLife
        @test Chemostats.get_curr_t(cell) ≈ log(2) atol = 1e-6
        @test Chemostats.get_curr_t(cell) ≈ Chemostats.get_curr_t(ref) atol = 1e-6
    end
end

@testset "DivideCallback aliasing" begin
    cb_eq = Chemostats.DivideCallback(jump_rn.V ~ jump_rn.V_d; interp_points=0, save_positions=(true, false))

    for _ in 1:3
        cell = run_to_division(hybrid_prob(cb_eq))
        @test Chemostats.get_state(cell) == Chemostats.CellState.EndOfLife
        @test Chemostats.get_curr_t(cell) ≈ log(2) atol = 1e-6
    end
end

@testset "DivideCallback with MTK equations (discrete=true)" begin
    Random.seed!(1)
    chem = Chemostat([ FlipCell() for _ in 1:20 ])
    Chemostats.simulate!(chem, 5.0, Chemostats.Direct())

    @test chem.status == Chemostats.ReturnCode.Success
    @test length(chem.pop) > 20   # population grew via repeated division
    @test all(c -> c.int[:X] <= 1, chem.pop)   # never overshoots past the divide point
end

@testset "SymbolicDivideContinuous / SymbolicDivideDiscrete" begin
    @variables X(t) = 0.0

    @test Chemostats.MTKDivideAffect() isa ModelingToolkit.ImperativeAffect

    sys_c = Chemostats.SymbolicDivideContinuous(X ~ 1.0; name = :divc)
    @test nameof(sys_c) == :divc
    cevs = ModelingToolkit.continuous_events(sys_c)
    @test length(cevs) == 1
    @test only(ModelingToolkit.equations(only(cevs))) == (X ~ 1.0)
    @test only(cevs).affect isa ModelingToolkit.ImperativeAffect

    sys_d = Chemostats.SymbolicDivideDiscrete(X >= 1; name = :divd)
    @test nameof(sys_d) == :divd
    devs = ModelingToolkit.discrete_events(sys_d)
    @test length(devs) == 1
    @test isequal(only(devs).conditions, X >= 1)
    @test only(devs).affect isa ModelingToolkit.ImperativeAffect
end
