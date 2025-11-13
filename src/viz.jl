using CairoMakie 

function plot(chem::Chemostat; tspan=(0, get_curr_t(chem)), Λ_gt=nothing, alg=nothing)
    fig = Figure()

    ax_N = Axis(fig[1,1], xlabel="Time", ylabel="N", yscale=log10)
    ax_Λ = Axis(fig[2,1], xlabel="Time", ylabel="Λ")
    ax_δ = Axis(fig[3,1], xlabel="Time", ylabel="δ")

    linkxaxes!([ax_N, ax_Λ, ax_δ])
    snaps = filter(snap -> tspan[1] <= snap.t <= tspan[2], chem.saved)
    tt = [ snap.t for snap in snaps ]
    NN = [ snap.N for snap in snaps ]

    lines!(ax_N, tt, NN)

    NN = [ est_N(snap) for snap in chem.saved ]
    lines!(ax_N, tt, NN)

    for f in [ 0., 0.1, 0.2 ]
        t_tgt = tspan[1] + (tspan[2] - tspan[1]) * f
        idx = findmin(t -> abs(t - t_tgt), tt)[2]
        t1 = tt[idx]
        ΛΛ = [ est_Λ(chem, t1, t) for t in tt ]
        lines!(ax_Λ, tt, ΛΛ)
    end

    isnothing(Λ_gt) || hlines!(ax_Λ, Λ_gt; linestyle=:dash, color=:black)

    if !isnothing(alg)
        δδ = [ get_δ(chem, t, alg) for t in tt ]
        lines!(ax_δ, tt, δδ)
        isnothing(Λ_gt) || hlines!(ax_δ, Λ_gt; linestyle=:dash, color=:black)
    end

    fig
end 