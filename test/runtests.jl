using Test
using Chemostats

# CI runs "Core" and "Extensions" as separate jobs; "Extensions" loads the heavy deps (MTK, Catalyst, Makie).
const GROUP = get(ENV, "GROUP", "All")

if GROUP in ("All", "Core")
	using Aqua

	@testset "Aqua.jl" begin
		Aqua.test_all(Chemostats)
	end

	@testset "Queue" begin include("queue.jl") end
	@testset "Deterministic" begin include("det.jl") end
	@testset "Unit" begin include("unit.jl") end
	@testset "Tree" begin include("tree.jl") end
	@testset "Scheduler" begin include("scheduler.jl") end
	@testset "Resilience" begin include("resilience.jl") end
	@testset "Growth rates" begin include("growth.jl") end
	@testset "Reinit" begin include("reinit.jl") end
end

if GROUP in ("All", "Extensions")
	@testset "MTK extension" begin include("mtk.jl") end
	@testset "Plotting extensions" begin include("plotting.jl") end
end

GROUP in ("All", "Core", "Extensions") || error("Unknown test GROUP: $GROUP")
