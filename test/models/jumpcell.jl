using JumpProcesses
using OrdinaryDiffEqTsit5
using SciMLBase
using SymbolicIndexingInterface: SymbolCache
using Chemostats

if !@isdefined(_MODELS_JUMPCELL_LOADED)
const _MODELS_JUMPCELL_LOADED = true

# Hybrid model: `P` is produced by the jump `0 --> P` at rate `k`, and `V` grows as dV/dt = V.
# u = [P, V], p = [k, V_d]; the `SymbolCache` lets tests index by name (`int[:P]`, `setp(int, :V_d)`).
const jumpcell_syms = SymbolCache([:P, :V], [:k, :V_d], :t)

function jumpcell_rhs!(du, u, p, t)
    du[1] = 0.0
    du[2] = u[2]
end

const jumpcell_birth = ConstantRateJump((u, p, t) -> p[1], int -> (int.u[1] += 1))

cb_cond_jump(u, t, int) = u[2] - int.p[2]
cb_div_jump = Chemostats.DivideCallback(cb_cond_jump; interp_points=0, save_positions=(true, false))

function divide_jumpcell(int)
    P = round(Int, int[:P])
    V_half = int[:V] / 2
    P1 = P ÷ 2
    P2 = P - P1

    [ (u0 = [:P => P1, :V => V_half], p = [:V_d => V_half + 1.0]),
      (u0 = [:P => P2, :V => V_half], p = [:V_d => V_half + 1.0]) ]
end

function JumpCell(alg = Tsit5(); k = 1.0, V_d = 2.0, remake = false)
    f = ODEFunction(jumpcell_rhs!; sys = jumpcell_syms)
    oprob = ODEProblem(f, [0.0, 1.0], (0., Inf), [k, V_d])
    prob = JumpProblem(oprob, Direct(), jumpcell_birth; callback = cb_div_jump)
    DECell(prob, alg, divide_jumpcell; remake)
end

jumpcell_aggregator(int) = int.opts.callback.discrete_callbacks[1].condition

end
