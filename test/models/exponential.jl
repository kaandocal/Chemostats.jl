using Random
using OrdinaryDiffEqTsit5
using Chemostats

f_exp(u, p, t) = -one(u)

cb_exp = Chemostats.DivideCallback((u, t, int) -> u; interp_points=0)

function divide_exp(int)
    p = int.p
    map(1:2) do _
        (u0 = randexp(), p = p)
    end
end

function ExponentialCell(p, alg = Tsit5())
    u0 = randexp()
    prob = ODEProblem(f_exp, u0, (0., 0.), p; callback = cb_exp)
    DECell(prob, alg, divide_exp)
end

cell_params(cell::DECell) = cell.int.p

function make_population(K; alg = Tsit5())
    [ ExponentialCell((; id), alg) for id in 1:K ]
end 

function make_stochastic_population(K; alg = Tsit5())
    [ StochasticCell((; id), alg) for id in 1:K ]
end 

function divide_stochastic(int)
    r = rand()
    p = int.p
    if r < 0.15
        nothing
    elseif r < 0.30
        [ (u0 = randexp(), p = p) ]
    else
        map(1:2) do _
            (u0 = randexp(), p = p)
        end
    end
end

function StochasticCell(p, alg = Tsit5())
    u0 = randexp()
    prob = ODEProblem(f_exp, u0, (0., 0.), p; callback = cb_exp)
    DECell(prob, alg, divide_stochastic)
end

nodivide(int) = nothing

function DyingCell(p, alg = Tsit5())
    u0 = randexp()
    prob = ODEProblem(f_exp, u0, (0., 0.), p; callback = cb_exp)
    DECell(prob, alg, nodivide)
end


make_dying_population(K; alg = Tsit5()) = [ DyingCell((; id), alg) for id in 1:K ]

struct CellError <: Exception
    id::Int
end
Base.showerror(io::IO, e::CellError) = print(io, "CellError: cell $(e.id) threw on purpose")

divide_error(int) = throw(CellError(int.p.id))

function ThrowingCell(p, alg = Tsit5())
    u0 = randexp()
    prob = ODEProblem(f_exp, u0, (0., 0.), p; callback = cb_exp)
    DECell(prob, alg, divide_error)
end

make_throwing_population(K; alg = Tsit5()) = [ ThrowingCell((; id), alg) for id in 1:K ]

has_cause(e, ::Type{T}) where T = e isa T
has_cause(e::CompositeException, ::Type{T}) where T = any(ex -> has_cause(ex, T), e.exceptions)
has_cause(e::TaskFailedException, ::Type{T}) where T = has_cause(e.task.result, T)
