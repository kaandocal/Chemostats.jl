module RecipesBaseExt

using Chemostats
using RecipesBase

@recipe function f(chem::Chemostats.Chemostat; Λ_gt=nothing, alg=nothing, tspan=nothing)
    tspan = isnothing(tspan) ? (0, chem.snaps[end].t) : tspan
    snaps = filter(snap -> tspan[1] <= snap.t <= tspan[2], chem.snaps)
    tt = [ snap.t for snap in snaps ]
    δ_panel = alg isa Chemostats.Lax

    layout := (δ_panel ? 4 : 3, 1)
    link := :x
    legend --> false

    @series begin
        subplot := 1
        ylabel := "log₁₀ N"
        tt, [ Chemostats.est_logN(snap) for snap in snaps ] ./ log(10)
    end

    @series begin
        subplot := 2
        ylabel := "simulated size L"
        tt, [ snap.N for snap in snaps ]
    end

    for f in [ 0., 0.1, 0.2 ]
        t_tgt = tspan[1] + (tspan[2] - tspan[1]) * f
        t1 = tt[findmin(t -> abs(t - t_tgt), tt)[2]]
        @series begin
            subplot := 3
            xlabel := δ_panel ? "" : "Time"
            ylabel := "growth rate Λ"
            tt, [ Chemostats.est_growth_rate(snaps, t1, t) for t in tt ]
        end
    end

    if !isnothing(Λ_gt)
        @series begin
            subplot := 3
            xlabel := δ_panel ? "" : "Time"
            ylabel := "growth rate Λ"
            seriestype := :path
            linestyle := :dash
            color := :black
            [ tt[1], tt[end] ], [ Λ_gt, Λ_gt ]
        end
    end

    if δ_panel
        @series begin
            subplot := 4
            xlabel := "Time"
            ylabel := "thinning rate δ"
            tt, [ Chemostats.get_thinning_rate(snaps, t, alg) for t in tt ]
        end
        if !isnothing(Λ_gt)
            @series begin
                subplot := 4
                xlabel := "Time"
                ylabel := "thinning rate δ"
                seriestype := :path
                linestyle := :dash
                color := :black
                [ tt[1], tt[end] ], [ Λ_gt, Λ_gt ]
            end
        end
    end
end

@recipe function f(snaps::AbstractVector{Chemostats.Snapshot})
    xlabel --> "Time"
    ylabel --> "log₁₀ N"
    legend --> false

    tt = [ snap.t for snap in snaps ]
    logNN = [ Chemostats.est_logN(snap) for snap in snaps ] ./ log(10)

    (tt, logNN)
end

end
