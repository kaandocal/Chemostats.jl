using Test 
using Chemostats

@testset "Deterministic (Serial)" begin global ensalg=Chemostats.EnsembleSerial(); include("det.jl") end
@testset "Deterministic (Threaded)" begin global ensalg=Chemostats.EnsembleThreads(); include("det.jl") end
#@testset "Yule" begin include("yule.jl") end
#@testset "Volume" begin include("vol.jl") end