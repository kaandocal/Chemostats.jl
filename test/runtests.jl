using Test 
using Chemostats
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
@testset "MTK extension" begin include("mtk.jl") end
