mutable struct DECell{I,DF,FT,P,Remake,R<:AbstractRNG}
    int::Union{I,Nothing}
    prob::P
    sol
    divide::DF
    state::CellState.T
    tshift::FT
    rng::R
end

"""
    DECell(prob, alg, divide; t0 = 0., reset_t=false, remake::Bool=false, rng=Xoshiro(rand(UInt64)))

Creates a new cell around a `SciMLBase.DEProblem`, simulated with algorithm `alg` starting at time `t0`.
The function `divide` takes a single integrator argument `int` and returns a list of offspring as a vector of `NamedTuples` with fields `u0` and `p`.

`reset_t` determines whether the integration time is started from `0` for every cell, which can avoid floating-point
errors for long simulations (many generations). If `reset_t` is set to `true`, each cell should keep track of its
starting time, e.g. via a parameter `t0` in `p`.

If `remake` is `true`, descendants will be created using `SciMLBase.remake`, which creates 
a copy of the underlying `AbstractProblem` (safest). Otherwise, descendants will be 
created using `SciMLBase.reinit`, which should be faster for `ODEProblem` and `SDEProblem` types.

**Note:** The problem should include a callback that calls `terminate!` when a cell reaches the end of its life.
"""
function DECell(prob, alg, divide; t0 = 0., reset_t=false, remake::Bool=false, rng=Xoshiro(rand(UInt64)))
    tshift, int = if reset_t
        int = copy_integrator(prob, alg, 0, rng)
        # Match daughters' tshift type (e.g. Int `t0` would otherwise change the cell type)
        convert(typeof(int.t), t0), int
    else
        nothing, copy_integrator(prob, alg, t0, rng)
    end

    DECell{typeof(int), typeof(divide), typeof(tshift), typeof(prob), remake, typeof(rng)}(
        int, prob, missing, divide, CellState.Newborn, tshift, rng
    )
end

get_cell_rng(cell::DECell) = cell.rng

# Independent stream of the same type for a new cell, derived deterministically from `rng`
spawn_rng(rng) = Random.seed!(copy(rng), rand(rng, UInt64))
# Not reproducible under multithreading; `copy` would turn it into a `Xoshiro`
spawn_rng(rng::TaskLocalRNG) = rng

# Callbacks can hold mutable state: give each cell fresh copies of `prob`'s callbacks, used exclusively
fresh_callback_kwargs(prob) = haskey(prob.kwargs, :callback) ?
    (; callback = deepcopy(prob.kwargs[:callback]), merge_callbacks = false) : (;)

# Copy to avoid cases where `init` aliases `p`, callbacks or a `JumpProblem`'s jump state
function copy_integrator(prob, alg, t0, rng)
    cb_kwargs = if prob isa SciMLBase.AbstractJumpProblem
        prob = deepcopy(prob)
        (;)
    else
        fresh_callback_kwargs(prob)
    end
    rng_kwargs = SciMLBase.supports_solve_rng(prob, alg) ? (; rng) : (;)

    int = init(prob, alg; tspan=(t0, Inf), cb_kwargs..., rng_kwargs...)
    int.p = deepcopy(int.p)
    int
end

process_int(int) = int.sol

function finalise!(cell::DECell)
    isnothing(cell.int) && return 

    cell.sol = process_int(cell.int)
    cell.int = nothing
end 

function reinit(old_int, prob; t0, rng, p=nothing, u0=nothing)
    int = copy_integrator(prob, old_int.alg, t0, rng)
    apply_reinit_overrides!(int, old_int; p, u0)
end

function reinit_remake(old_int, prob; t0, rng, p=nothing, u0=nothing)
    base_prob = SciMLBase.remake(prob; u0 = copy(old_int.u), p = old_int.p)
    int = copy_integrator(base_prob, old_int.alg, t0, rng)
    apply_reinit_overrides!(int, old_int; p, u0)
end

function apply_reinit_overrides!(int, old_int; p, u0)
    if !isnothing(p)
        if p isa AbstractVector{<:Pair}
            setp = SymbolicIndexingInterface.setp(int, first.(p))
            setp(int, last.(p))
        else
            int.p = deepcopy(p)
        end
    end

    u = copy(old_int.u)

    if !isnothing(u0) && !isempty(u0)
        if first(u0) isa Pair
            setu = SymbolicIndexingInterface.setu(int, first.(u0))
            setu(u, last.(u0))
        else
            u = copy(u0)
        end
    end

    SciMLBase.reinit!(int, u; reinit_dae=false)

    int
end

# Dispatches on a cell's `remake` type parameter
reinit_dispatch(remake::Bool, old_int, prob; t0, rng, p=nothing, u0=nothing) =
    remake ? reinit_remake(old_int, prob; t0, rng, p, u0) : reinit(old_int, prob; t0, rng, p, u0)

function get_children(parent::DECell{I,DF,FT,P,Remake,R}, p=nothing) where {I,DF,FT,P,Remake,R}
    @assert !isnothing(parent.int)

    children = parent.divide(parent.int)
    isnothing(children) && return typeof(parent)[]

    t = get_curr_t(parent)

    map(children) do ch
        tshift, t0 = if isnothing(parent.tshift)
            nothing, t
        else
            t, zero(t)
        end

        rng = spawn_rng(parent.rng)
        int = reinit_dispatch(Remake, parent.int, parent.prob; t0, rng, u0=ch.u0, p=ch.p)
        DECell{typeof(int), DF, typeof(tshift), P, Remake, R}(
            int, parent.prob, missing, parent.divide, CellState.Newborn, tshift, rng
        )
    end
end

function get_curr_t(cell::DECell) 
    ret = if !isnothing(cell.int)
        cell.int.t 
    else 
        last(cell.sol.t)
    end 

    if !isnothing(cell.tshift)
        ret += cell.tshift 
    end 

    ret 
end 

get_state(cell::DECell) = cell.state

function set_state!(cell::DECell, state::CellState.T)
    savevalues!(cell)
    cell.state = state
end 

is_alive(cell::DECell) = get_state(cell) == CellState.Newborn || get_state(cell) == CellState.Alive

function init_cell!(cell::DECell)
    @argcheck cell.state == CellState.Newborn "Cannot reinitialise cell in state $(cell.state)"
    
    set_state!(cell, CellState.Alive)
end

function kill!(cell::DECell, t=get_curr_t(cell))
    @assert !isnothing(cell.int)

    cell.state = CellState.Killed

    SciMLBase.done(cell.int) || terminate!(cell.int)
    finalise!(cell)
end

function die!(cell::DECell)
    @assert !isnothing(cell.int)

    if get_state(cell) != CellState.EndOfLife
        @warn "Called `die!` on cell in state $(get_state(cell))"
        return 
    end 

    set_state!(cell, CellState.Dead)
    finalise!(cell)
end 

function divide!(cell::DECell)
    @assert !isnothing(cell.int)

    if get_state(cell) != CellState.EndOfLife
        @warn "Cell in state $(get_state(cell)) tried to divide"
        return 
    end 

    set_state!(cell, CellState.Divided)
    finalise!(cell)
end

savevalues!(cell::DECell) = savevalues!(cell.int)

function step!(cell::DECell, dt, p)
    @assert !isnothing(cell.int)

    isnothing(p) || @warn "Parameter arguments to DECell are ignored" maxlog=1

    if get_state(cell) != CellState.Alive
        @warn "Tried to simulate cell in state $(get_state(cell))"
        return 
    end 

    step!(cell.int, dt, true)
    savevalues!(cell)

    if SciMLBase.done(cell.int)
        if !SciMLBase.successful_retcode(cell.int.sol)
            @warn "Cell solver errored, removing cell..."
            set_state!(cell, CellState.Killed)
            finalise!(cell)
        else 
            set_state!(cell, CellState.EndOfLife)
        end
    end 
end

function clone_cell(cell::DECell{I,DF,FT,P,Remake,R}, t) where {I,DF,FT,P,Remake,R}
    @assert !isnothing(cell.int)

    tshift = isnothing(cell.tshift) ? zero(t) : cell.tshift
    @check cell.int.sol.t[1] <= t - tshift <= get_curr_t(cell)

    rng = spawn_rng(cell.rng)
    int = reinit_dispatch(Remake, cell.int, cell.prob; t0=t - tshift, rng, u0=cell.int.sol(t - tshift), p=cell.int.p)
    DECell{typeof(int), DF, FT, P, Remake, R}(int, cell.prob, missing, cell.divide, CellState.Newborn, cell.tshift, rng)
end

###

"""
    DivideCallback(condition; discrete=false, kwargs...)

Implements a differential equation callback to check whether a cell has reached the end of its lifetime.
If `discrete=false` (default), returns a `ContinuousCallback` that calls `terminate!` once `condition` 
equals zero. If `discrete=true`, returns a `DiscreteCallback`.  
"""
function DivideCallback(condition; discrete::Bool=false, kwargs...)
    if discrete
        SciMLBase.DiscreteCallback(condition, SciMLBase.terminate!; kwargs...)
    else
        SciMLBase.ContinuousCallback(condition, SciMLBase.terminate!; kwargs...)
    end
end
