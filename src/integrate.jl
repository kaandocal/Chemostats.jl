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

function PopIntegrator(chem::Chemostat, alg::AbstractAlgorithm, ensalg::EnsembleAlgorithm; tstops = Float64[])
    t0 = get_curr_t(chem)

    tstops = filter(t -> t0 < t, tstops)
    tstops = unique(sort(tstops))

    queue = create_queue(chem.pop, ensalg) 
    
    PopIntegrator(chem, alg, queue, Float64.(tstops), t0, t0, t0, 
                  chem.saved[end].nsim, chem.saved[end].log_f, 
                  SciMLBase.ReturnCode.Default)
end 

Snapshot(int::PopIntegrator) = Snapshot(int.t, length(int.queue), int.nsim, int.log_f)
savevalues!(int::PopIntegrator) = push!(int.chem.saved, Snapshot(int))

function add_tstop!(int::PopIntegrator, t)
    @argcheck t >= int.t_next "Cannot add tstop at time $t: middle of simulation"

    t in int.tstops && return 

    push!(int.tstops, t)
    sort!(int.tstops)
end 

function simulate!(chem::Chemostat, tmax, alg::AbstractAlgorithm, 
                   ensalg::EnsembleAlgorithm = EnsembleSerial(); saveat=[ tmax ], kwargs...)
    int = PopIntegrator(chem, alg, ensalg; tstops=saveat)
    init!(int.alg, int)
    simulate!(int, tmax, ensalg)
    chem
end


function simulate!(int::PopIntegrator, tmax, ensalg::EnsembleAlgorithm; kwargs...)
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

    while true
        if isempty(queue) || int.retcode != ReturnCode.Default 
            break
        elseif length(queue) > Nmax  
            @warn "Population size exceeds $Nmax, aborting. Consider adjusting Nmax."
            int.retcode = ReturnCode.MaxIters
            break 
        end 

        cell = first(queue)
        get_curr_t(cell) >= int.t_next && break

        pop!(queue)
        if get_state(cell) == CellState.Newborn 
            init_cell!(cell)
            int.nsim += 1
        end 

        if get_state(cell) == CellState.Alive 
            process_cell!(int, cell, int.t_next; δ, kwargs...)
        elseif get_state(cell) == CellState.Divided
            process_division!(int, cell; kwargs...)
        end

        update_queue!(int, int.alg, get_curr_t(cell))
    end

    int.log_f += δ * (int.t_next - t0)

    int.t = if int.retcode == ReturnCode.MaxIters 
        first(queue).t
    else 
        int.t_next 
    end 

    save && savevalues!(int)

    int
end



function step!(int::PopIntegrator, tmax, ::EnsembleThreads; Nmax=Int(1e7), save=false, kwargs...)
    @unpack chem, queue = int 
    out = create_queue(empty(chem.pop), EnsembleThreads())
    
    int.t_next = find_next_t(int, tmax)
    t0 = int.t 
    int.t_next <= t0 && return int

    δ = get_δ(int, int.alg)
    @sync for i in 1:Threads.nthreads()
        Threads.@spawn begin
            register_listener!(queue)
            # skip this 
            # update_queue!(int, int.alg)

            while true 
                if length(queue) > Nmax  
                    @warn "Population size exceeds $Nmax, aborting. Consider adjusting Nmax."
                    lock(int.queue) do 
                        int.retcode = ReturnCode.MaxIters
                    end
                    break 
                end 

                cell = fetch!(queue)
                if isnothing(cell) || int.retcode != ReturnCode.Default 
                    break
                elseif get_curr_t(cell) >= int.t_next
                    #@info "Pushing cell at time $(get_curr_t(cell)) out..."
                    push!(out, cell)
                    continue
                end 

                if get_state(cell) == CellState.Newborn 
                    #@info "At time $(get_curr_t(cell)): Newborn cell"
                    init_cell!(cell)
                    int.nsim += 1
                end 

                if get_state(cell) == CellState.Alive 
                    #@info "At time $(get_curr_t(cell)): Alive cell"
                    process_cell!(int, cell, int.t_next; δ, kwargs...)
                elseif get_state(cell) == CellState.Divided
                    #@info "At time $(get_curr_t(cell)): Dividing cell"
                    process_division!(int, cell; kwargs...)
                end

                #@info "End at time $(get_curr_t(cell))"
            end
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
    step!(cell, tmax - tb, int.chem.model, int.chem.env)
    t = get_curr_t(cell)

    @check t <= tmax + 1e-6
    @check get_state(cell) != CellState.Alive || t > tb "Cell simulation did not increase time"

    # Cells are culled with rate δ
    if δ > 0 
        t_kill = tb + randexp() / δ

        if t_kill < t
            kill_cell!(cell, t_kill)
            save_leaves && push!(int.chem.leaves, cell)
            return
        end
    end

    # This uses a threadsafe call
    push!(int.queue, cell)
end 

function process_division!(int::PopIntegrator, cell; save_lineages=false, kwargs...)
    @argcheck get_state(cell) == CellState.Divided 
    #@info "Dividing cell..."

    # Cell dies or divides
    offspr = get_offspring(int, cell; save_lineages)
    offspr_filtered = filter_offspring(offspr, int.alg)

    lockqueue(int.queue) do 
        N_full = length(int.queue) + length(offspr)
        N_obs = length(int.queue) + length(offspr_filtered)

        N_full / N_obs
        _append!(int.queue, offspr_filtered)
        int.log_f += log(N_full) - log(N_obs)
    end
end

function resize_pop!(int, L::Int, t)
    lockqueue(int.queue) do 
        if isempty(int.queue)
            @warn "No cells left in chemostat, terminating..."
            int.retcode = ReturnCode.Unstable
            return
        end

        N_start = length(int.queue)
        
        while length(int.queue) < L
            _duplicate_random!(int.queue, t)
        end

        _truncate_queue!(int.queue, L)

        int.log_f += log(N_start) - log(L)
    end
end