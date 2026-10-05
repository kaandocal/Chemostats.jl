using Catalyst
using ModelingToolkit
using JumpProcesses
using OrdinaryDiffEqTsit5
using SciMLBase
using Chemostats

jump_rn = @reaction_network begin
    @species P(t) = 0
    @variables V(t) = 1.

    @parameters k V_d

    @equations begin
        D(V) ~ V
    end

    k, 0 --> P
end

_jumpcell_ref_prob = HybridProblem(jump_rn, [:P => 0, :V => 1.0], (0., Inf), [:k => 1.0, :V_d => 2.0])
getV_jump = ModelingToolkit.getu(_jumpcell_ref_prob, :V)

cb_cond_jump(u, t, int) = getV_jump(u) - int.ps[:V_d]
cb_div_jump = Chemostats.DivideCallback(cb_cond_jump; interp_points=0, save_positions=(true, false))

function divide_jumpcell(int)
    P = round(Int, int[:P])
    V_half = int[:V] / 2
    P1 = P ÷ 2
    P2 = P - P1

    [ (u0 = [:P => P1, :V => V_half], p = [:V_d => V_half + 1.0]),
      (u0 = [:P => P2, :V => V_half], p = [:V_d => V_half + 1.0]) ]
end

function JumpCell(alg = Tsit5(); k = 1.0, V_d = 2.0)
    prob = HybridProblem(jump_rn, [:P => 0, :V => 1.0], (0., Inf), [:k => k, :V_d => V_d]; callback = cb_div_jump)
    DECell(prob, alg, divide_jumpcell)
end

jumpcell_aggregator(int) = int.opts.callback.discrete_callbacks[1].condition
