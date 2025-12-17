struct Snapshot 
    t::Float64
    N::Int 
    nsim::Int           # number of (partially) simulated cells
    log_f::Float64      # (log) fraction of cells retained
end 

est_logN(snap::Snapshot) = log(snap.N) + snap.log_f
est_N = exp ∘ est_logN

function est_Λ(before::Snapshot, after::Snapshot) 
    before.t >= after.t && return NaN
    (est_logN(after) - est_logN(before)) / (after.t - before.t)
end

function snapshot_str(snap::Snapshot)
    @unpack t, N, nsim = snap
    N_est = round(est_N(snap); sigdigits=3)

    return "t=$t, N=$N/$N_est, nsim=$nsim"
end 

function Base.show(io::IO, snap::Snapshot)
    stats = snapshot_str(snap)
    print(io, "Snapshot($stats)")
end

### 

struct Chemostat{C,P}
    pop::Vector{C}
    p::P
    saved::Vector{Snapshot}
    leaves::Vector{C}
    all_cells::Vector{C}               # all cells 
end

function Chemostat(cells, p=nothing)
    saved = [ Snapshot(0., length(cells), 0, 0.) ]
    Chemostat(deepcopy(cells), p, saved, empty(cells), empty(cells))
end

function Base.show(io::IO, chem::Chemostat)
    stats = snapshot_str(chem.saved[end])
    print(io, "Chemostat($stats)")
end

const SnapshotVec = Union{AbstractVector{Snapshot},Chemostat}

get_curr_t(chem::Chemostat) = get_curr_t(chem.saved)
get_curr_t(snaps::AbstractVector{Snapshot}) = snaps[end].t 

get_snapshot(chem::Chemostat, args...) = get_snapshot(chem.saved, args...)

function get_snapshot(snaps::AbstractVector{Snapshot}, t = snaps[end].t)
    @assert issorted(snaps; by = snap -> snap.t)
    i = searchsorted(map(snap -> snap.t, snaps), t)
    isempty(i) && throw(ArgumentError("No snapshot saved at time t=$t"))
    snaps[last(i)]
end 

""" estimate the real population size by extrapolation """
est_logN(snaps::SnapshotVec, t = get_curr_t(snaps)) = est_logN(get_snapshot(snaps, t))

""" estimate the growth rat between time t0 and t """
est_Λ(snaps::SnapshotVec, t0, t) = est_Λ(get_snapshot(snaps, t0), get_snapshot(snaps, t))
est_Λ(snaps::SnapshotVec, t = get_curr_t(snaps)) = est_Λ(snaps, zero(t), t)

# # Population history stuff 

# add_lineage!(set, ::Missing) = set
# function add_lineage!(set, cell::Cell)
#     add_lineage!(set, cell.anc)
#     union!(set, [cell])
# end 

# function collect_cells!(chem::Chemostat)
#     for cell in chem.pop
#         add_lineage!(chem.all_cells, cell)
#     end 

#     for cell in chem.leaves 
#         add_lineage!(chem.all_cells, cell)
#     end
# end

# filter_cells(cells, t) = filter(cell -> cell.t[1] <= t < cell.t[end], cells)

# function all_cells(chem::Chemostat)
#     isempty(chem.all_cells) && collect_cells!(chem)
#     chem.all_cells 
# end

# get_all_cells(chem::Chemostat, t::Float64) = filter_cells(get_all_cells(chem), t)
