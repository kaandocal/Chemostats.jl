mutable struct PopIntegrator{CT <: Chemostat, A <: AbstractAlgorithm, QT <: QueueType}
    chem::CT
    alg::A
    queue::QT
    tstops::Vector{Float64}
    t0::Float64
    t::Float64
    t_next::Float64
    nsim::Int
    log_f::Float64
    retcode::ReturnCode.T
end

function default_alg(alg::AbstractAlgorithm) 
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

    queue = create_queue(chem.pop, ensalg) 
    
    PopIntegrator(chem, alg, queue, Float64.(tstops), t0, t0, t0, 
                  chem.snaps[end].nsim, chem.snaps[end].log_f, 
                  SciMLBase.ReturnCode.Default)
end 

Snapshot(int::PopIntegrator) = Snapshot(int.t, length(int.queue), int.nsim, int.log_f)
savevalues!(int::PopIntegrator) = push!(int.chem.snaps, Snapshot(int))

function add_tstop!(int::PopIntegrator, t)
    @argcheck t >= int.t_next "Cannot add tstop at time $t: middle of simulation"

    t in int.tstops && return 

    push!(int.tstops, t)
    sort!(int.tstops)
end 

""" 
    simulate!(chem, tmax, alg, ensalg = EnsembleThreads(); saveat = [])

Simulates the chemostat until time `tmax` with algorithm `alg` (see [Simulation Algorithms](@ref) for a list of algorithms).
`ensalg` can be `EnsembleSerial()` for single-threaded simulations, or `EnsembleThreads()` for multithreading (if `alg` support multithreading).
The default `ensalg` is `EnsembleThreads()` for multithreaded algorithms, or `EnsembleSerial()` for single-threaded algorithms.
"""    
function simulate!(chem::Chemostat, tmax, alg::AbstractAlgorithm, 
                   ensalg::EnsembleAlgorithm = default_alg(alg); saveat=[ tmax ], kwargs...)
    int = PopIntegrator(chem, alg, ensalg; tstops=saveat)
    init!(int.alg, int)
    simulate!(int, tmax, ensalg; kwargs...)
    chem
end


function simulate!(int::PopIntegrator, tmax, ensalg::EnsembleAlgorithm = default_alg(alg); kwargs...)
    while int.t < tmax && int.retcode == ReturnCode.Default
        update_algorithm!(int.alg, int)
        int.retcode == ReturnCode.Default || break
        step!(int, tmax, ensalg; save=true, kwargs...)
    end

    empty!(int.chem.pop)
    extract_queue!(int.chem.pop, int.queue)

    if int.retcode == ReturnCode.Default 
        int.retcode = ReturnCode.Success
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

function worker_task(int::PopIntegrator, out::QueueType; Nmax=Int(1e7), δ=0., kwargs...)
    register_listener!(int.queue)

    while true 
        if length(int.queue) > Nmax  
            @warn "Population size exceeds $Nmax, aborting. Consider adjusting Nmax."
            lock(int.queue) do 
                int.retcode = ReturnCode.MaxIters
            end
            break 
        end 

        cell = fetch!(int.queue)
        if isnothing(cell) || int.retcode != ReturnCode.Default 
            break
        elseif get_curr_t(cell) >= int.t_next
            push!(out, cell)
            continue
        end 

        if get_state(cell) == CellState.Newborn 
            init_cell!(cell)
            int.nsim += 1
        end 

        if get_state(cell) == CellState.Alive 
            process_cell!(int, cell, int.t_next; δ, kwargs...)
        elseif get_state(cell) == CellState.EndOfLife
            process_eol!(int, cell; kwargs...)
        end

        # We assume this is threadsafe (`Strict` does not support multithreading)
        update_queue!(int, int.alg, get_curr_t(cell))
    end
end 

function step!(int::PopIntegrator, tmax, ensalg::Union{EnsembleSerial,EnsembleThreads}; save=false, kwargs...)
    @unpack chem, queue = int 
    out = create_queue(empty(chem.pop), ensalg)
    
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

    int.t = if int.retcode == ReturnCode.MaxIters 
        first(queue).t
    else 
        int.t_next 
    end 

    int.queue = out
    save && savevalues!(int)

    int
end

function process_cell!(int::PopIntegrator, cell, tmax; δ=0., save_leaves=false, kwargs...)
    @argcheck get_state(cell) == CellState.Alive 

    tb = get_curr_t(cell)
    step!(cell, tmax - tb, int.chem.p)
    t = get_curr_t(cell)

    @check get_state(cell) != CellState.Alive || t > tb "Cell simulation did not increase time"

    # Cells are culled with rate δ
    if δ > 0 
        t_kill = tb + randexp() / δ

        if t_kill < t
            kill!(cell, t_kill)
            save_leaves && push!(int.chem.leaves, cell)
            return
        end
    end

    # This uses a threadsafe call
    push!(int.queue, cell)
end 

function process_eol!(int::PopIntegrator, cell; save_lineages=false, kwargs...)
    @argcheck get_state(cell) == CellState.EndOfLife
    @debug "Dividing cell..."

    # Cell dies or divides
    offspr = get_offspring(int, cell; save_lineages)

    if isempty(offspr)
        die!(cell)
        return 
    end 

    offspr_filtered, Δlog_f = filter_offspring(offspr, int.alg)
    divide!(cell)

    lock(int.queue) do 
        @debug "Appending $(length(offspr_filtered))/$(length(offspr)) cells..."
        _append!(int.queue, offspr_filtered)
        int.log_f += Δlog_f
    end
end

function _resize_pop!(int, L::Int, t)
    if isempty(int.queue)
        @warn "No cells left in chemostat, terminating..."
        int.retcode = ReturnCode.Unstable
        return
    end

    N_start = length(int.queue)
    
    while length(int.queue) < L
        _clone_random!(int.queue, t)
    end

    _truncate_queue!(int.queue, L)

    int.log_f += log(N_start) - log(L)
end
