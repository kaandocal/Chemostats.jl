using SafeTestsets

# CI runs "Trees", "Core" and "Extensions" in turn; "Extensions" loads the heavy deps (MTK, Catalyst, Makie).
const GROUP = get(ENV, "GROUP", "All")
GROUP in ("All", "Trees", "Core", "Extensions") || error("Unknown test GROUP: $GROUP")

if GROUP in ("All", "Trees")
    @time @safetestset "Trees" include("Trees/runtests.jl")
end

if GROUP in ("All", "Core")
    @time @safetestset "Aqua.jl" begin
        using Aqua, Chemostats
        Aqua.test_all(Chemostats)
    end

    @time @safetestset "Queue" include("queue.jl")
    @time @safetestset "Chemostat" include("chemostat.jl")
    @time @safetestset "Algorithms" include("algorithms.jl")
    @time @safetestset "DECell" include("decell.jl")
    @time @safetestset "Offspring" include("offspring.jl")
    @time @safetestset "Recording" include("recording.jl")
    @time @safetestset "Resilience" include("resilience.jl")
    @time @safetestset "Growth rates" include("growth.jl")
    @time @safetestset "RNG" include("rng.jl")
end

if GROUP in ("All", "Extensions")
    @time @safetestset "MTK extension" include("mtk.jl")
    @time @safetestset "Plotting extensions" include("plotting.jl")
end
