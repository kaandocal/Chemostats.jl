using Test
using Chemostats

# CI runs "Trees", "Core" and "Extensions" in turn; "Extensions" loads the heavy deps (MTK, Catalyst, Makie).
const GROUP = get(ENV, "GROUP", "All")
GROUP in ("All", "Trees", "Core", "Extensions") || error("Unknown test GROUP: $GROUP")

if GROUP in ("All", "Trees")
	# Own module, so the Trees test helpers don't leak into the other groups
	@eval module TreesTests include(joinpath(@__DIR__, "Trees", "runtests.jl")) end
end

if GROUP in ("All", "Core")
	using Aqua

	@testset "Aqua.jl" begin
		Aqua.test_all(Chemostats)
	end

	@testset "Queue" begin include("queue.jl") end
	@testset "Chemostat" begin include("chemostat.jl") end
	@testset "Algorithms" begin include("algorithms.jl") end
	@testset "DECell" begin include("decell.jl") end
	@testset "Offspring" begin include("offspring.jl") end
	@testset "Recording" begin include("recording.jl") end
	@testset "Resilience" begin include("resilience.jl") end
	@testset "Growth rates" begin include("growth.jl") end
	@testset "RNG" begin include("rng.jl") end
end

if GROUP in ("All", "Extensions")
	@testset "MTK extension" begin include("mtk.jl") end
	@testset "Plotting extensions" begin include("plotting.jl") end
end
