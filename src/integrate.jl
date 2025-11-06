const TOrder = Base.By(int -> int.t)

mutable struct PopIntegrator{CT <: Chemostat, CI <: CellIntegrator, A <: AbstractAlgorithm}
    chem::CT
    alg::A
    queue::BinaryHeap{CI}
    tstops::Vector{Float64}
    t::Float64
    t_next::Float64
    nsim::Int
    log_f::Float64
    retcode::ReturnCode.T
end

function PopIntegrator(chem::Chemostat, alg::AbstractAlgorithm; tstops = Float64[])
    t0 = get_t(chem)

    tstops = filter(t -> t0 < t, tstops)
    tstops = unique(sort(tstops))

    queue = BinaryHeap(TOrder, [ CellIntegrator(cell) for cell in chem.pop ])
    
    PopIntegrator(chem, alg, queue, Float64.(tstops), t0, t0, 
                  chem.saved[end].nsim, chem.saved[end].log_f, 
                  SciMLBase.ReturnCode.Default)
end 

Snapshot(int::PopIntegrator) = Snapshot(int.t, length(int.queue), int.nsim, int.log_f)

function add_tstop!(int::PopIntegrator, t)
    @argcheck t >= int.t_next "Cannot add tstop before current t"

    t in int.tstops && return 

    push!(int.tstops, t)
    sort!(int.tstops)
end 

function init!(int::PopIntegrator)
    init!(int.alg, int)
    push!(int.chem.saved, Snapshot(int))
end 

function simulate!(chem::Chemostat, tmax, alg::AbstractAlgorithm; saveat=[ tmax ], kwargs...)
    int = PopIntegrator(chem, alg; tstops=saveat)
    init!(int)
    solve!(int, tmax)
    chem
end

function solve!(int::PopIntegrator, tmax; kwargs...)
    isnothing(int.chem.env) || simulate_env!(int.chem.env, tmax)

    while int.t < tmax && int.retcode == ReturnCode.Default
        update_algorithm!(int.alg, int)
        int.retcode == ReturnCode.Default || break
        step!(int, tmax; save=true, kwargs...)
    end

    empty!(int.chem.pop)
    while !isempty(int.queue)
        cell = pop!(int.queue)
        #@check cell.sol.status == CellState.Alive
        push!(int.chem.pop, cell.sol)
    end

    if int.retcode == ReturnCode.Default 
        int.retcode = ReturnCode.Success
    end 
end 

function find_next_t(int::PopIntegrator, tmax)
    idx = findfirst(t -> t > int.t, int.tstops)
    idx == nothing && return tmax 
    min(int.tstops[idx], tmax)
end 

function step!(int::PopIntegrator, tmax; Nmax=Int(1e7), save=false, kwargs...)
    @unpack chem, queue = int 
    
    int.t_next = find_next_t(int, tmax)
    t0 = int.t 
    int.t_next <= t0 && return int

    δ = get_δ(int.alg)

    while !isempty(queue)
        if length(queue) > Nmax  
            @warn "Population size exceeds $Nmax, aborting. Consider adjusting Nmax."
            int.retcode = ReturnCode.MaxIters
            break 
        end 

        cell = first(queue)
        cell.t >= int.t_next && break

        pop!(queue)
        if cell.sol.status == CellState.Alive 
            process_cell!(int, cell, int.t_next; δ, kwargs...)
        elseif cell.sol.status == CellState.Divided
            process_division!(int, cell; kwargs...)
        else 
            process_death!(int, cell; kwargs...)
        end
    end

    int.log_f += δ * (int.t_next - t0)

    int.t = if int.retcode == ReturnCode.MaxIters 
        first(queue).t
    else 
        int.t_next 
    end 

    save && push!(int.chem.saved, Snapshot(int))

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

    offspr = [ CellIntegrator(Cell(anc, cell.t, u, p)) for (u, p) in offspr ]
    N_expected = length(int.queue) + length(offspr) 
    int.nsim += add_offspring!(int.queue, offspr, int.alg)
    int.log_f += log(N_expected) - log(length(int.queue))
end

function process_death!(int::PopIntegrator, cell; save_lineages=false, kwargs...)
    @argcheck cell.sol.status == CellState.Dead || cell.sol.status == CellState.Killed
    
    N_expected = length(int.queue)
    add_offspring!(int.queue, [], int.alg)
    int.log_f += log(N_expected) - log(length(int.queue))
end