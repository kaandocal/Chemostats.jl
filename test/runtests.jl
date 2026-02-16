using Test 
using Chemostats
using Aqua

@testset "Aqua.jl" begin
	Aqua.test_all(Chemostats)
end

@testset "Deterministic" begin include("det.jl") end
@testset "Growth rates" begin include("growth.jl") end
