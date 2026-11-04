# Observables are user-defined functions `f(node, a)` of a node at age `a` (time since birth). Per-cell
# summaries are functions `f(node)` of a [`Node`](@ref), and population
# statistics are (weighted) means over node collections (a tree, `view(tree, ids)` or `ancestors(tree, id)`).

# Integrates `f(node, age)` over ages in `[a, b]` using adaptive quadrature
quadgk_integral(f, node, a, b) = first(quadgk(age -> f(node, age), a, b))

"""
    time_average(f, node::Node; tspan = (-Inf, Inf), integrate)
    time_average(f; tspan = (-Inf, Inf), integrate)

Returns the time average of `f(node, age)` over the cell's lifetime, truncated to lie within `tspan`. Returns `nothing` for nodes without data, or for cells not alive within `tspan`.

`time_average(f; kwargs...)` is a convenience wrapper for `node -> time_average(f, node; kwargs...)`.

Time averages are computed by integration, using the callback `integrate(f, node, a, b)`. By default, this uses adaptive Gaussian quadrature ([QuadGK.jl](@extref QuadGK :doc:`index`)).
"""
function time_average(f, node::Node; tspan = (-Inf, Inf), integrate=quadgk_integral)
    isnothing(node.data) && return nothing
    a = max(node.t_birth, tspan[1]) - node.t_birth
    b = min(node.t_end, tspan[2]) - node.t_birth
    b > a ? integrate(f, node, a, b) / (b - a) : nothing
end

time_average(f; kwargs...) = node -> time_average(f, node; kwargs...)

###

"""
    time_weight(node::Node; tspan = (-Inf, Inf))
    time_weight(; tspan = (-Inf, Inf))

Weights cells by their lifetime within `tspan` (living cells up to their last update). Used to compute time averages in
lineages and populations. `time_weight(; tspan)` returns the function `node -> time_weight(node; tspan)`.
"""
time_weight(node::Node; tspan = (-Inf, Inf)) = max(0., min(node.t_end, tspan[2]) - max(node.t_birth, tspan[1]))
time_weight(; tspan = (-Inf, Inf)) = node -> time_weight(node; tspan)

# Weighted mean and weighted sum of squared deviations (`nothing` without values), total weight and number of values
# of `f(node)` over the nodes of `nodes` with non-zero weight, skipping nodes where `f` returns `nothing` (West's algorithm)
function weighted_moments(f, nodes, w)
    mS, W, n = nothing, 0., 0
    for node in values(nodes)
        w_i = w(node)
        iszero(w_i) && continue

        v = f(node)
        isnothing(v) && continue

        n += 1
        W += w_i
        if isnothing(mS)
            mS = (float(v), zero(float(v)))
        else
            m, S = mS
            δ = v .- m
            m = m .+ (w_i / W) .* δ
            mS = (m, S .+ w_i .* δ .* (v .- m))
        end
    end

    mS, W, n
end

"""
    mean(f, nodes)
    mean(f, nodes, w)

Returns the mean of `f(node)` over all nodes in `nodes`, optionally weighted by `w`. The default returns the average over cells; if `w = time_weight`, returns time averages. Ignores nodes where `f` returns `nothing`.

For example:

```julia
# Time average of an observable `f(node, age)` over the population within `tspan`
mean(time_average(f; tspan), tree, time_weight(; tspan))

# Mean lifetime of the cells along the backward lineage of node `id`
mean(lifetime, ancestors(tree, id))
```
"""
function StatsBase.mean(f, nodes::NodeCollection, w = Returns(1.))
    mS, _, _ = weighted_moments(f, nodes, w)
    isnothing(mS) ? NaN : first(mS)
end

"""
    var(f, nodes; corrected=true)
    var(f, nodes, w; corrected=true)

Returns the variance of `f(node)` over all nodes in `nodes`, optionally weighted by `w`, see [`mean`](@ref). If `corrected=true`, applies Bessel's correction (for probability weights).

**Note:** This is the variance of per-cell summaries. The variance of per-cell averages is smaller than the variance across the population (law of total variance).
"""
function StatsBase.var(f, nodes::NodeCollection, w = Returns(1.); corrected::Bool=true)
    last(mean_and_var(f, nodes, w; corrected))
end

"""
    mean_and_var(f, nodes; corrected=true)
    mean_and_var(f, nodes, w; corrected=true)

Returns the mean and variance of `f(node)` over all nodes in `nodes` in a single pass, optionally weighted by `w`, see [`mean`](@ref) and [`var`](@ref). If `corrected=true`, applies Bessel's correction to the variance (for probability weights).
"""
function StatsBase.mean_and_var(f, nodes::NodeCollection, w = Returns(1.); corrected::Bool=true)
    mS, W, n = weighted_moments(f, nodes, w)
    isnothing(mS) && return NaN, NaN
    m, S = mS
    n <= corrected && return m, NaN
    m, corrected ? S ./ W .* (n / (n - 1)) : S ./ W
end

###

"""
    snapshot(f, nodes, t)

Returns `f(node, age)` for all nodes in `nodes` alive at time `t`. Nodes without data are skipped. Requires a sorted tree.
"""
function snapshot(f, nodes, t)
    [ f(node, t - node.t_birth) for node in alive_at(nodes, t) if !isnothing(node.data) ]
end

"""
    snapshot(f, nodes, ts::AbstractVector)

Returns the concatenated snapshots of the nodes in `nodes` at the times `ts`. Requires a sorted tree.
"""
snapshot(f, nodes, ts::AbstractVector) = collect(Iterators.flatten(snapshot(f, nodes, t) for t in ts))

###

# Records: age-conditioned and per-cycle values over node collections. Nodes without data (e.g. living cells) are skipped.

"""
    at_age(f, nodes, a)

Returns `f(node, a)` for the nodes in `nodes` that have reached age `a`. Nodes without data, or nodes that do not reach age `a`, are skipped.
"""
function at_age(f, nodes, a)
    [ f(node, a) for node in values(nodes) if !isnothing(node.data) && lifetime(node) >= a ]
end

"""
    at_birth(f, nodes)

Returns `f(node, 0)` for the nodes of the node collection `nodes`.
"""
at_birth(f, nodes) = at_age(f, nodes, 0)

"""
    at_division(f, nodes)

Returns `f(node, τ)` for the nodes of the node collection `nodes` that divided, where `τ` is their generation time.
Requires a sorted tree.
"""
function at_division(f, nodes)
    ensure_sorted(nodes)
    tree = get_tree(nodes)
    [ f(node, lifetime(node)) for (id, node) in pairs(nodes) if !isnothing(node.data) && has_divided(tree, id) ]
end
