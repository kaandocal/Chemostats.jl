using Test 

@testset "Deterministic" begin include("det.jl") end
@testset "Yule" begin include("yule.jl") end
@testset "Volume" begin include("vol.jl") end