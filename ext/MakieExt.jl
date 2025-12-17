module MakieExt

using Chemostats, Makie

Makie.plot(chem::Chemostat; kwargs...) = Makie.plot(chem.saved; kwargs...)

function Makie.plot(snaps::AbstractVector{Chemostats.Snapshot};
                    tspan=(0, snaps[end].t), Λ_gt=nothing, alg=nothing)
    fig = Figure()

    tt = [ snap.t for snap in snaps ]

    ax_N = Axis(fig[1,1], xlabel="Time", ylabel="log10(N)")
    ax_Λ = Axis(fig[2,1], xlabel="Time", ylabel="Λ")
    ax_δ = if alg isa Chemostats.Lax
        Axis(fig[3,1], xlabel="Time", ylabel="δ")
    else 
        nothing 
    end

    linkxaxes!([ax_N, ax_Λ])
    snaps = filter(snap -> tspan[1] <= snap.t <= tspan[2], snaps)

    logNN = [ Chemostats.est_logN(snap) for snap in snaps ] ./ log(10)
    lines!(ax_N, tt, logNN)

    for f in [ 0., 0.1, 0.2 ]
        t_tgt = tspan[1] + (tspan[2] - tspan[1]) * f
        idx = findmin(t -> abs(t - t_tgt), tt)[2]
        t1 = tt[idx]
        ΛΛ = [ Chemostats.est_Λ(snaps, t1, t) for t in tt ]
        lines!(ax_Λ, tt, ΛΛ)
    end

    isnothing(Λ_gt) || hlines!(ax_Λ, Λ_gt; linestyle=:dash, color=:black)

    if alg isa Chemostats.Lax
        linkxaxes!([ax_N, ax_δ])
        δδ = [ Chemostats.get_δ(snaps, t, alg) for t in tt ]
        lines!(ax_δ, tt, δδ)
        isnothing(Λ_gt) || hlines!(ax_δ, Λ_gt; linestyle=:dash, color=:black)
    end

    fig
end 

end
