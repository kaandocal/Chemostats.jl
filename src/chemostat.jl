struct Snapshot 
    t::Float64
    N::Int 
    nsim::Int           # number of (partially) simulated cells
    log_f::Float64      # (log) fraction of cells retained
end 

est_logN(snap::Snapshot) = log(snap.N) + snap.log_f
est_N = exp ∘ est_logN

struct Chemostat{C,M,E}
    pop::Vector{C}
    model::M
    env::E
    saved::Vector{Snapshot}
    leaves::Vector{C}
    all_cells::Vector{C}               # all cells 
end

function Base.show(io::IO, chem::Chemostat)
    compact = get(io, :compact, true)
    t = chem.saved[end].t
    N = chem.saved[end].N
    N_est = est_N(chem)
    if compact 
        print(io, "Chemostat(t=$t, N=$N/$N_est)")
    else 
        println(io, "Chemostat:")
        println(io, " t = $t")
        print(io, " Cells = $N (out of ≈ $N_est)")
    end 
end

function Chemostat(cells, model, env)
    saved = [ Snapshot(0., length(cells), length(cells), 0.) ]
    Chemostat(deepcopy(cells), model, env, saved, empty(cells), empty(cells))
end

get_t(chem::Chemostat) = chem.saved[end].t 

""" estimate the real population size by extrapolation """

function est_logN(chem::Chemostat, t = get_t(chem))
    i = findfirst(snap -> snap.t == t, chem.saved)
    if isnothing(i) 
        throw(ArgumentError("No snapshot saved at time t=$t"))
    end

    est_logN(chem.saved[i])
end

function est_Λ(chem::Chemostat, t0, t) 
    t0 >= t && return NaN
    (est_logN(chem, t) - est_logN(chem, t0)) / (t - t0)
end

est_Λ(chem::Chemostat, t = get_t(chem)) = est_Λ(chem, zero(t), t)

# Population history stuff 

add_lineage!(set, ::Missing) = set
function add_lineage!(set, cell::Cell)
    add_lineage!(set, cell.anc)
    union!(set, [cell])
end 

function collect_cells!(chem::Chemostat)
    for cell in chem.pop
        add_lineage!(chem.all_cells, cell)
    end 

    for cell in chem.leaves 
        add_lineage!(chem.all_cells, cell)
    end
end

# filter_cells(cells, t) = filter(cell -> cell.t[1] <= t < cell.t[end], cells)

# function all_cells(chem::Chemostat)
#     isempty(chem.all_cells) && collect_cells!(chem)
#     chem.all_cells 
# end

# get_all_cells(chem::Chemostat, t::Float64) = filter_cells(get_all_cells(chem), t)
