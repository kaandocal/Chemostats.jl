using Test 
using Chemostats

@testset "Deterministic (Serial)" let ensalg=Chemostats.EnsembleSerial(); include("det.jl") end
@testset "Deterministic (Threaded)" let ensalg=Chemostats.EnsembleThreaded(); include("det.jl") end
#@testset "Yule" begin include("yule.jl") end
#@testset "Volume" begin include("vol.jl") end