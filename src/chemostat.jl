"""
    struct Snapshot
        t::Float64
        N::Int 
        nsim::Int
        log_f::Float64
    end 

Saves the state of a simulate dpopulation at a fixed time `t`. Here `N` is the observed size of the population,
and `log_f` is minus the log fraction of the true population observed. This allows us to use subsampling algorithms
(see [Simulation Algorithms](@ref)) that only track a small subset of a full population and extrapolate.

The true population size can be estimated using [`est_logN`](@ref) or [`est_N`](@ref). The field `nsim` 
counts the number of cells simulated until the current snapshot (including partially simulated and removed cells).
"""
struct Snapshot 
    t::Float64
    N::Int 
    nsim::Int           # number of (partially) simulated cells
    log_f::Float64      # (log) fraction of cells retained
end 

"""
    est_logN(snap::Snapshot)
    
Estimates the true log-population size from a snapshot.
"""
est_logN(snap::Snapshot) = log(snap.N) + snap.log_f

"""
    est_N(snap::Snapshot)
    
Estimates the true population size from a snapshot.
"""
est_N = exp ∘ est_logN

"""
    est_Λ(snap1::Snapshot, snap2::Snapshot)

Estimates the growth rates from two snapshots at times ``t_1`` and ``t_2`` via 

``\\hat \\Lambda(t_1, t_2) = \\frac{\\log \\hat N(t_2) - \\log \\hat N(t_1)}{t_2 - t_1}``

Here ``\\hat N(t_1)`` and ``\\hat N(t_2)`` are estimated via [`est_N`](@ref).
"""
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

"""
    struct Chemostat{C,P}
        pop::Vector{C}
        p::P
        snaps::Vector{Snapshot}
    end    

Chemostat object. The field `pop` contains the current population of cells.
The field `p` can be a user-defined parameter object (defaults to `nothing`).
`snaps` contains a list of `Snapshot`s saved at different times.
"""
struct Chemostat{C,P}
    pop::Vector{C}
    p::P
    snaps::Vector{Snapshot}
end

function Chemostat(cells, p=nothing; copy=false)
    snaps = [ Snapshot(0., length(cells), 0, 0.) ]
    Chemostat(copy ? deepcopy(cells) : cells, p, snaps)
end

function Base.show(io::IO, chem::Chemostat)
    stats = snapshot_str(chem.snaps[end])
    print(io, "Chemostat($stats)")
end

const SnapshotVec = Union{AbstractVector{Snapshot},Chemostat}

get_curr_t(chem::Chemostat) = get_curr_t(chem.snaps)
get_curr_t(snaps::AbstractVector{Snapshot}) = snaps[end].t 

get_snapshot(chem::Chemostat, args...) = get_snapshot(chem.snaps, args...)

"""
    get_snapshot(chem::Chemostat, t = chem.t)
    get_snapshot(snaps::AbstractVector{Snapshot}, t = snaps[end].t)

Searches and returns for a population snapshot saved at time `t`. Errors if no snapshots are found.
If multiple snapshots at time `t` exist, this function returns the last one. 
"""
function get_snapshot(snaps::AbstractVector{Snapshot}, t = snaps[end].t; atol=1e-6)
    @assert issorted(snaps; by = snap -> snap.t)
    i = searchsortedlast(map(snap -> snap.t, snaps), t + atol)
    if i < 1 || snaps[i].t < t - atol
        throw(ArgumentError("No snapshot saved at time t=$t"))
    end 
    
    snaps[last(i)]
end 


"""
    est_logN(snaps::AbstractVector{Snapshot}, t)

Estimates the log-size of a population at time `t` by combining [`get_snapshot`](@ref) and [`est_logN`](@ref).
`snaps` can alternatively be a `Chemostat` instance.
"""
est_logN(snaps::SnapshotVec, t = get_curr_t(snaps)) = est_logN(get_snapshot(snaps, t))

"""
    est_Λ(snaps::AbstractVector{Snapshot}, t1 = 0., t2)

Estimates the growth rate of a population between times t1 and t2 
by combining [`get_snapshot`](@ref) and [`est_Λ`](@ref). `snaps` can alternatively be a `Chemostat` instance.
"""
est_Λ(snaps::SnapshotVec, t0, t) = est_Λ(get_snapshot(snaps, t0), get_snapshot(snaps, t))
est_Λ(snaps::SnapshotVec, t = get_curr_t(snaps)) = est_Λ(snaps, zero(t), t)