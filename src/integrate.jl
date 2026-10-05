mutable struct PopIntegrator{CT <: Chemostat, A <: AbstractAlgorithm, QT <: ThreadedQueue}
    chem::CT
    alg::A
    queue::QT
    tstops::Vector{Float64}
    t0::Float64
    t::Float64
    t_next::Float64
    nsim::Threads.Atomic{Int}
    log_f::Float64
    @atomic retcode::ReturnCode.T
    tree_lock::ReentrantLock
end

function default_ensalg(alg::AbstractAlgorithm) 
    if is_parallel(alg) && Threads.nthreads() > 1
        EnsembleThreads()
    else 
        EnsembleSerial()
    end 
end

function PopIntegrator(chem::Chemostat, alg::AbstractAlgorithm, ensalg::EnsembleAlgorithm; tstops = Float64[])
    t0 = get_curr_t(chem)

    tstops = filter(t -> t0 < t, tstops)
    tstops = unique(sort(tstops))

    queue = ThreadedQueue(chem.pop) 
    
    PopIntegrator(chem, alg, queue, Float64.(tstops), t0, t0, t0,
                  Threads.Atomic{Int}(chem.snaps[end].nsim), chem.snaps[end].log_f,
                  SciMLBase.ReturnCode.Default, ReentrantLock())
end 

Snapshot(int::PopIntegrator) = Snapshot(int.t, length(int.queue), int.nsim[], int.log_f)
savevalues!(int::PopIntegrator) = push!(int.chem.snaps, Snapshot(int))

function add_tstop!(int::PopIntegrator, t)
    @argcheck t >= int.t_next "Cannot add tstop at time $t: middle of simulation"

    t in int.tstops && return 

    push!(int.tstops, t)
    sort!(int.tstops)
end 

"""
    simulate!(chem, tmax, alg, ensalg = EnsembleThreads(); saveat = [], Nmax = 1e7, throw_on_error = false)

Simulates the chemostat until time `tmax` with algorithm `alg` (see [Simulation Algorithms](@ref) for a list of algorithms).
`ensalg` can be `EnsembleSerial()` for single-threaded simulations, or `EnsembleThreads()` for multithreading (if `alg` support multithreading).
The default `ensalg` is `EnsembleThreads()` for multithreaded algorithms, or `EnsembleSerial()` for single-threaded algorithms.

If the population size exceeds `Nmax`, the simulation exits (useful for `Direct` and `Lax`).

If `throw_on_error` is `true`, exceptions that occur while simulating a cell will stop the simulation. 
If it is `false`, the exception will be caught and the cell will be removed from the population.
"""
function simulate!(chem::Chemostat, tmax, alg::AbstractAlgorithm,
                   ensalg::EnsembleAlgorithm = default_ensalg(alg); saveat=[ tmax ], throw_on_error::Bool=false, kwargs...)
    int = PopIntegrator(chem, alg, ensalg; tstops=saveat)
    init!(int.alg, int)
    simulate!(int, tmax, ensalg; throw_on_error, kwargs...)
    chem
end

function simulate!(int::PopIntegrator, tmax, ensalg::EnsembleAlgorithm = default_ensalg(int.alg); throw_on_error::Bool=false, kwargs...)
    while int.t < tmax && int.retcode == ReturnCode.Default
        update_algorithm!(int.alg, int)
        int.retcode == ReturnCode.Default || break
        step!(int, tmax, ensalg; save=true, throw_on_error, kwargs...)

        simplify!(int.chem.tree, iter_unsafe(int.queue))
    end

    empty!(int.chem.pop)
    extract_queue!(int.chem.pop, int.queue)

    if int.retcode == ReturnCode.Default
        @atomic int.retcode = ReturnCode.Success
    end

    int
end 

### TODO: sorted search
### This must be strict (>) for Lax to work
function find_next_t(int::PopIntegrator, t=int.t)
    idx = findfirst(s -> s > t, int.tstops)
    idx == nothing && return Inf
    int.tstops[idx]
end 

function step!(int::PopIntegrator, tmax, ensalg::EnsembleAlgorithm; kwargs...)
    error("Ensemble algorithm $ensalg not supported")
end 

function worker_task(int::PopIntegrator, out::ThreadedQueue; Nmax=Int(1e7), δ=0., throw_on_error::Bool=false, kwargs...)
    register_listener!(int.queue)

    while true
        try
            if length(int.queue) > Nmax
                @warn "Population size exceeds $Nmax, aborting. Consider adjusting Nmax."

                @atomicreplace int.retcode ReturnCode.Default => ReturnCode.MaxIters
                release_and_notify!(int.queue)
                return
            end

            cell = fetch!(int.queue)
            isnothing(cell) && return

            if int.retcode != ReturnCode.Default
                push!(int.queue, cell)
                release_and_notify!(int.queue)
                return
            elseif get_curr_t(cell) >= int.t_next
                push!(out, cell)
                continue
            end

            if get_state(cell) == CellState.Newborn
                init_cell!(cell)
                Threads.atomic_add!(int.nsim, 1)
            end

            if get_state(cell) == CellState.Alive
                process_cell!(int, cell, int.t_next; δ, kwargs...)
            elseif get_state(cell) == CellState.EndOfLife
                process_eol!(int, cell; kwargs...)
            end

            # We assume this is threadsafe (`Strict` does not support multithreading)
            update_queue!(int, int.alg, get_curr_t(cell))
        catch e
            # Simulating the cell caused an error
            if throw_on_error
                @atomicreplace int.retcode ReturnCode.Default => ReturnCode.Failure
                release_and_notify!(int.queue)
                rethrow()
            else
                showerror(stderr, e, catch_backtrace())
                flush(stderr)
            end
        end
    end
end

function step!(int::PopIntegrator, tmax, ensalg::Union{EnsembleSerial,EnsembleThreads}; save=false, kwargs...)
    @unpack chem, queue = int 
    out = ThreadedQueue(empty(chem.pop))
    
    if ensalg isa EnsembleThreads && !is_parallel(int.alg)
        error("Algorithm $(typeof(int.alg)) does not support parallelisation")
    end 

    δ = get_δ(int, int.alg)
    int.t_next = min(find_next_t(int), tmax)
    t0 = int.t 
    int.t_next <= t0 && return int

    if ensalg isa EnsembleSerial
        worker_task(int, out; δ, kwargs...)
    elseif ensalg isa EnsembleThreads
        @sync for i in 1:Threads.nthreads()
            Threads.@spawn worker_task(int, out; δ, kwargs...)
        end
    end

    int.log_f += δ * (int.t_next - t0)

    if int.retcode == ReturnCode.Default
        int.t = int.t_next
        int.queue = out
    else
        int.t = isempty(queue) ? int.t_next : get_curr_t(first(queue))

        merged = empty(chem.pop)
        extract_queue!(merged, queue)
        extract_queue!(merged, out)
        int.queue = ThreadedQueue(merged)
    end

    save && savevalues!(int)

    int
end

function process_cell!(int::PopIntegrator, cell, tmax; δ=0., kwargs...)
    @argcheck get_state(cell) == CellState.Alive 

    tb = get_curr_t(cell)
    step!(cell, tmax - tb, int.chem.p)
    t = get_curr_t(cell)

    @check get_state(cell) != CellState.Alive || t > tb "Cell simulation did not increase time"

    if get_state(cell) == CellState.Killed 
        @lock int.tree_lock add_leaf!(int.chem.tree, cell)
        return
    end 

    # Cells are culled with rate δ
    if δ > 0 
        t_kill = tb + randexp() / δ

        if t_kill < t
            kill!(cell, t_kill)
            @lock int.tree_lock add_leaf!(int.chem.tree, cell)
            return
        end
    end

    # This uses a threadsafe call
    push!(int.queue, cell)
end 

function process_eol!(int::PopIntegrator, cell; kwargs...)
    @argcheck get_state(cell) == CellState.EndOfLife
    @debug "Dividing cell..."

    # Cell dies or divides
    children = get_children(cell, int.chem.p)

    if isempty(children)
        die!(cell)
        @lock int.tree_lock add_leaf!(int.chem.tree, cell)
        return 
    end 

    children_filtered, Δlog_f = filter_offspring(children, int.alg)
    divide!(cell)
    @lock int.tree_lock add_offspring!(int.chem.tree, cell, Tuple(children_filtered))

    lock(int.queue) do 
        @debug "Appending $(length(children_filtered))/$(length(children)) cells..."
        _append!(int.queue, children_filtered)
        int.log_f += Δlog_f
    end
end

function _resize_pop_unsafe!(int, L::Int, t)
    if isempty(int.queue)
        @warn "No cells left in chemostat, terminating..."
        @atomic int.retcode = ReturnCode.Unstable
        return
    end

    N_start = length(int.queue)

    while length(int.queue) > L 
        j = rand(1:length(int.queue))
        cell = _popat!(int.queue, j)
        kill!(cell, t)
        @lock int.tree_lock add_leaf!(int.chem.tree, cell)
    end 
    
    while length(int.queue) < L
        _clone_random!(int.queue, t)
    end

    int.log_f += log(N_start) - log(L)
end

