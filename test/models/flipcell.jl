using Catalyst
using ModelingToolkit
using JumpProcesses
using OrdinaryDiffEqTsit5
using Chemostats

flip_rn = @reaction_network begin
    @species X(t) = 0
    @variables V(t) = 1.
    @parameters k
    @equations begin
        D(V) ~ 0
    end
    k, 0 --> X
end

cb_flip = Chemostats.DivideCallback(flip_rn.X ~ 1.0; discrete=true)

function divide_flipcell(int)
    [ (u0 = [:X => 0, :V => 1.0], p = [:k => int.ps[:k]]) for _ in 1:2 ]
end

function FlipCell(; k = 1.0)
    prob = HybridProblem(flip_rn, [:X => 0, :V => 1.0], (0., Inf), [:k => k]; callback = cb_flip)
    DECell(prob, Tsit5(), divide_flipcell)
end
