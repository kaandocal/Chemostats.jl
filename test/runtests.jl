using Test 
using Chemostats

@testset "Deterministic (Serial)" begin ensalg=Chemostats.EnsembleSerial(); include("det.jl") end
@testset "Deterministic (Threaded)" begin ensalg=Chemostats.EnsembleThreads(); include("det.jl") end
#@testset "Yule" begin include("yule.jl") end
#@testset "Volume" begin include("vol.jl") end