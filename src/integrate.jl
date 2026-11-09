mutable struct PopIntegrator{CT <: Chemostat, A <: AbstractAlgorithm, QT <: ThreadedQueue, R <: AbstractRNG}
    chem::CT
    alg::A
    queue::QT
    rng::R                      # population-level random choices (e.g. by `Strict`)
    tstops::Vector{Float64}
    saveat::Vector{Float64}     # populations are marked in the tree at these times
    t0::Float64
    t::Float64
    t_next::Float64
    @atomic nsim::Int
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

function PopIntegrator(chem::Chemostat, alg::AbstractAlgorithm, ensalg::EnsembleAlgorithm; tstops = Float64[],
                       rng::AbstractRNG = Xoshiro(rand(UInt64)))
    t0 = get_curr_t(chem)

    tstops = filter(t -> t0 < t, tstops)
    tstops = unique(sort(tstops))

    queue = ThreadedQueue(chem.pop) 
    
    PopIntegrator(chem, alg, queue, rng, Float64.(tstops), Float64.(tstops), t0, t0, t0,
                  chem.snaps[end].nsim, chem.snaps[end].log_f,
                  SciMLBase.ReturnCode.Default, ReentrantLock())
end 

Snapshot(int::PopIntegrator) = Snapshot(int.t, length(int.queue), int.nsim, int.log_f)
savevalues!(int::PopIntegrator) = push!(int.chem.snaps, Snapshot(int))

# Moves the cells in `out` back into `int.queue`, e.g. after early termination
function merge_queues!(int::PopIntegrator, out::ThreadedQueue)
    merged = empty(int.chem.pop)
    extract_queue!(merged, int.queue)
    extract_queue!(merged, out)
    int.queue = ThreadedQueue(merged)
end

function add_tstop!(int::PopIntegrator, t)
    @argcheck t >= int.t_next "Cannot add tstop at time $t: middle of simulation"

    t in int.tstops && return 

    push!(int.tstops, t)
    sort!(int.tstops)
end 

"""
    simulate(cells, tmax, alg, ensalg = EnsembleThreads(); kwargs...)

Utility function that creates a [`Chemostat`](@ref) with initial population `cells`
and runs [`simulate!`](@ref). Returns the created `Chemostat` object.
"""
function simulate(cells, args...; p=nothing, copy=false, save_mode::Symbol=:none, kwargs...)
    chem = Chemostat(cells; p, copy, save_mode)
    simulate!(chem, args...; kwargs...)
end 

"""
    simulate!(chem, tmax, alg, ensalg = EnsembleThreads(); saveat = [], Nmax = 1e7, throw_on_error = false)

Simulates the chemostat until time `tmax` with algorithm `alg` (see [Simulation Algorithms](@ref) for a list of algorithms).
`ensalg` can be `EnsembleSerial()` for single-threaded simulations, or `EnsembleThreads()` for multithreading (if `alg` support multithreading).
The default `ensalg` is `EnsembleThreads()` for multithreaded algorithms, or `EnsembleSerial()` for single-threaded algorithms.

If the population size exceeds `Nmax`, the simulation exits (useful for `Direct` and `Lax`).

If a population tree is recorded, cells alive at times `saveat` are saved together with their 
ancestors.

If `throw_on_error` is `true`, exceptions that occur while simulating a cell will stop the simulation. 
If it is `false`, the exception will be caught and the cell will be removed from the population.
"""
function simulate!(chem::Chemostat, tmax, alg::AbstractAlgorithm,
                   ensalg::EnsembleAlgorithm = default_ensalg(alg); saveat=[ tmax ], throw_on_error::Bool=false,
                   rng::AbstractRNG = Xoshiro(rand(UInt64)), kwargs...)
    int = PopIntegrator(chem, alg, ensalg; tstops=saveat, rng)
    try
        init!(int.alg, int)
        int.t0 in saveat && mark_tree!(chem)
        simulate!(int, tmax, ensalg; throw_on_error, kwargs...)
    catch
        # Keep the cells simulated so far
        empty!(chem.pop)
        extract_queue!(chem.pop, int.queue)
        chem.status = ReturnCode.Failure
        rethrow()
    finally
        # Also after early termination or an exception
        update_tree!(chem.tree)
    end
    chem.status = int.retcode
    chem
end

function simulate!(int::PopIntegrator, tmax, ensalg::EnsembleAlgorithm = default_ensalg(int.alg); throw_on_error::Bool=false, kwargs...)
    while int.t < tmax && int.retcode == ReturnCode.Default
        update_algorithm!(int.alg, int)
        int.retcode == ReturnCode.Default || break
        step!(int, tmax, ensalg; save=true, throw_on_error, kwargs...)

        prune_tree!(int.chem)
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
    throw(ArgumentError("Ensemble algorithm $ensalg not supported, use EnsembleSerial() or EnsembleThreads()"))
end 

function worker_task(int::PopIntegrator, cell, out::ThreadedQueue; Nmax=Int(1e7), δ=0., throw_on_error::Bool=false, kwargs...)
    t_start = get_curr_t(cell)
    processed = false
    try
        if int.retcode != ReturnCode.Default
            push!(int.queue, cell)
            return false
        end

        if length(int.queue) > Nmax
            @warn "Population size exceeds $Nmax, aborting. Consider adjusting Nmax."
            @atomicreplace int.retcode ReturnCode.Default => ReturnCode.MaxIters
            push!(int.queue, cell)
            return false
        end

        if get_curr_t(cell) >= int.t_next
            push!(out, cell)
            return
        end

        if get_state(cell) == CellState.Newborn
            init_cell!(cell)
            @atomic int.nsim += 1
        end

        if get_state(cell) == CellState.Alive
            process_cell!(int, cell, int.t_next; δ, kwargs...)
        elseif get_state(cell) == CellState.EndOfLife
            process_eol!(int, cell; kwargs...)
        end
        processed = true

        # We assume this is threadsafe (`Strict` does not support multithreading)
        update_queue!(int, int.alg, get_curr_t(cell))
    catch e
        if e isa CellException
            # The cell is removed from the population
            kill!(cell)
            @lock int.tree_lock record_drop!(int.chem.tree, cell)
        end

        if throw_on_error || !(e isa CellException)
            # Return the popped cell to the population, unless it was passed on or removed (e.g. killed above)
            !processed && get_state(cell) in (CellState.Newborn, CellState.Alive, CellState.EndOfLife) &&
                push!(int.queue, cell)

            @atomicreplace int.retcode ReturnCode.Default => ReturnCode.Failure
            rethrow()
        else
            showerror(stderr, e, catch_backtrace())
            flush(stderr)
            # Replace the cell when it was taken from the queue, before any other cell moved on
            update_queue!(int, int.alg, t_start)
        end
    end
end

function step!(int::PopIntegrator, tmax, ensalg::Union{EnsembleSerial,EnsembleThreads}; save=false, kwargs...)
    (; chem, queue) = int
    out = ThreadedQueue(empty(chem.pop))

    if ensalg isa EnsembleThreads && !is_parallel(int.alg)
        throw(ArgumentError("Algorithm $(typeof(int.alg)) does not support parallelisation, use EnsembleSerial()"))
    end

    δ = get_δ(int, int.alg)
    int.t_next = min(find_next_t(int), tmax)
    t0 = int.t
    int.t_next <= t0 && return int

    wp = WorkerPool(queue)

    try
        if ensalg isa EnsembleSerial
            run!(wp) do cell
                worker_task(int, cell, out; δ, kwargs...)
            end
        elseif ensalg isa EnsembleThreads
            @sync for i in 1:Threads.nthreads()
                Threads.@spawn run!(wp) do cell
                    worker_task(int, cell, out; δ, kwargs...)
                end
            end
        end
    catch
        # Cells in `out` would otherwise be lost
        merge_queues!(int, out)
        rethrow()
    end

    if int.retcode == ReturnCode.Default
        int.t = int.t_next
        int.queue = out
        int.t in int.saveat && mark_tree!(chem)
    else
        int.t = isempty(queue) ? int.t_next : get_curr_t(first(queue))
        merge_queues!(int, out)
    end

    # Only up to the time reached, which is earlier on early termination
    int.log_f += δ * (int.t - t0)

    save && savevalues!(int)

    int
end

# Interrupts and setup errors stop the simulation instead of removing the cell
stops_simulation(e) = e isa InterruptException || e isa NonShareableError

function process_cell!(int::PopIntegrator, cell, tmax; δ=0., kwargs...)
    @argcheck get_state(cell) == CellState.Alive 

    tb = get_curr_t(cell)
    try
        step!(cell, tmax - tb, int.chem.p)
    catch e
        stops_simulation(e) && rethrow()
        throw(CellException(e))
    end
    t = get_curr_t(cell)

    @check get_state(cell) != CellState.Alive || t > tb "Cell simulation did not increase time"

    if get_state(cell) == CellState.Killed 
        @lock int.tree_lock record_end!(int.chem.tree, cell, t, get_data(cell))
        return
    end 

    # Cells are culled with rate δ
    if δ > 0 
        t_kill = tb + randexp(get_cell_rng(cell)) / δ

        if t_kill < t
            data = get_data(cell)
            kill!(cell, t_kill)
            @lock int.tree_lock record_end!(int.chem.tree, cell, t_kill, data)
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
    children = try
        get_children(cell, int.chem.p)
    catch e
        stops_simulation(e) && rethrow()
        throw(CellException(e))
    end

    t, data = get_curr_t(cell), get_data(cell)

    if isempty(children)
        die!(cell)
        @lock int.tree_lock record_end!(int.chem.tree, cell, t, data)
        return 
    end 

    children_filtered, Δlog_f = filter_offspring(children, int.alg, get_cell_rng(cell))
    divide!(cell)
    # Before the children are queued, so their own events come after their birth
    @lock int.tree_lock record_division!(int.chem.tree, cell, children_filtered, t, data)

    lock(int.queue) do 
        @debug "Appending $(length(children_filtered))/$(length(children)) cells..."
        _append!(int.queue, children_filtered)
        int.log_f += Δlog_f
    end
end
