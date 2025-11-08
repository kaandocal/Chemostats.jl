const TOrder = Base.By(int -> int.t)

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

create_queue(cells, ::EnsembleSerial) = BinaryHeap(TOrder, [ CellIntegrator(cell) for cell in cells ])
create_queue(cells, ::EnsembleThreads) = ThreadedBinaryHeap(create_queue(cells, EnsembleSerial()))

function extract_queue!(pop, queue::BinaryHeap)
    while !isempty(queue)
        cell = pop!(queue)
        push!(pop, cell.sol)
    end
end 

function extract_queue!(pop, queue::ThreadedBinaryHeap)
    lockqueue!(queue)
    extract_queue!(pop, queue.heap)
    unlockqueue!(queue)
end 

function PopIntegrator(chem::Chemostat, alg::AbstractAlgorithm, ensalg::EnsembleAlgorithm; tstops = Float64[])
    t0 = get_t(chem)

    tstops = filter(t -> t0 < t, tstops)
    tstops = unique(sort(tstops))

    queue = create_queue(chem.pop, ensalg) 
    
    PopIntegrator(chem, alg, queue, Float64.(tstops), t0, t0, t0, 
                  chem.saved[end].nsim, chem.saved[end].log_f, 
                  SciMLBase.ReturnCode.Default)
end 

simulate_env!(::Nothing, tmax) = nothing 

savevalues!(int::PopIntegrator) = push!(int.chem.saved, Snapshot(int))

Snapshot(int::PopIntegrator) = Snapshot(int.t, length(int.queue), int.nsim, int.log_f)

function add_tstop!(int::PopIntegrator, t)
    @argcheck t >= int.t_next "Cannot add tstop at time $t: middle of simulation"

    t in int.tstops && return 

    push!(int.tstops, t)
    sort!(int.tstops)
end 

function init!(int::PopIntegrator)
    init!(int.alg, int)
    savevalues!(int)
end 

function simulate!(chem::Chemostat, tmax, alg::AbstractAlgorithm, 
                   ensalg::EnsembleAlgorithm = EnsembleSerial(); saveat=[ tmax ], kwargs...)
    int = PopIntegrator(chem, alg, ensalg; tstops=saveat)
    init!(int)
    solve!(int, tmax, ensalg)
    chem
end

function solve!(int::PopIntegrator, tmax, ensalg::EnsembleAlgorithm; kwargs...)
    simulate_env!(int.chem.env, tmax)

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
end 

function find_next_t(int::PopIntegrator, tmax)
    idx = findfirst(t -> t > int.t, int.tstops)
    idx == nothing && return tmax 
    min(int.tstops[idx], tmax)
end 

function step!(int::PopIntegrator, tmax, ensalg::EnsembleAlgorithm; kwargs...)
    error("Ensemble algorithm $ensalg not supported")
end 

function step!(int::PopIntegrator, tmax, ::EnsembleSerial; Nmax=Int(1e7), save=false, kwargs...)
    @unpack chem, queue = int 
    
    int.t_next = find_next_t(int, tmax)
    t0 = int.t 
    int.t_next <= t0 && return int

    δ = get_δ(int, int.alg)

    # THREAD SAFE
    while true
        should_sync(queue, int.alg) && sync!(int, int.alg)

        if isempty(queue) || int.retcode != ReturnCode.Default 
            break
        elseif length(queue) > Nmax  
            @warn "Population size exceeds $Nmax, aborting. Consider adjusting Nmax."
            int.retcode = ReturnCode.MaxIters
            break 
        end 

        cell = first(queue)
        cell.t >= int.t_next && break

        pop!(queue)
        if cell.sol.status == CellState.Newborn 
            init_cell!(cell)
            int.nsim += 1
        end 

        if cell.sol.status == CellState.Alive 
            process_cell!(int, cell, int.t_next; δ, kwargs...)
        elseif cell.sol.status == CellState.Divided
            process_division!(int, cell; kwargs...)
        end
    end
    # END THREAD

    int.log_f += δ * (int.t_next - t0)

    int.t = if int.retcode == ReturnCode.MaxIters 
        first(queue).t
    else 
        int.t_next 
    end 

    save && savevalues!(int)

    int
end

function process_cell!(int::PopIntegrator, cell, tmax; δ=0., save_leaves=false, kwargs...)
    @argcheck cell.sol.status == CellState.Alive 

    tb = cell.t 
    step!(cell, tmax - tb, int.chem.model, int.chem.env)

    @check cell.t <= tmax + 1e-9
    @check cell.sol.status != CellState.Alive || cell.t > tb "Cell simulation for time 0"

    # Cells are culled with rate δ
    if δ > 0 
        t_kill = tb + rand(Exponential(1 / δ))

        if t_kill < cell.t
            kill_cell!(cell, t_kill)
            save_leaves && push!(int.chem.leaves, cell.sol)
            return
        end
    end

    # This is threadsafe
    push!(int.queue, cell)
end 

function process_division!(int::PopIntegrator, cell; save_lineages=false, kwargs...)
    @argcheck cell.sol.status == CellState.Divided 

    # Cell dies or divides
    offspr = get_offspring(cell, int.chem.model, int.chem.env)

    anc = if save_lineages 
        cell.sol 
    else 
        missing 
    end 

    offspr = map(((u, p),) -> CellIntegrator(Cell(anc, cell.t, u, p)), offspr)

    N_full = length(int.queue) + length(offspr)
    offspr_filtered = filter_offspring(offspr, int.alg)
    N_obs = length(int.queue) + length(offspr_filtered)

    lockqueue(int.queue) do 
        append!(int.queue, offspr_filtered)
        int.log_f += log(N_full) - log(N_obs)
    end
end