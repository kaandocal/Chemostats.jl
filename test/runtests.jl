using Test 
using Chemostats
using Aqua

@testset "Aqua.jl" begin
	Aqua.test_all(Chemostats)
end

#@testset "Deterministic (Serial)" begin global ensalg=Chemostats.EnsembleSerial(); include("deterministic.jl") end
#@testset "Deterministic (Threads)" begin global ensalg=Chemostats.EnsembleThreads(); include("deterministic.jl") end

#@testset "Yule (Serial)" begin global ensalg=Chemostats.EnsembleSerial(); include("yule.jl") end
#@testset "Yule (Threads)" begin global ensalg=Chemostats.EnsembleThreads(); include("yule.jl") end

@testset "Two-type Markov (Serial)" begin global ensalg=Chemostats.EnsembleSerial(); include("twotype.jl") end
@testset "Two-type Markov (Threads)" begin global ensalg=Chemostats.EnsembleThreads(); include("twotype.jl") end

#@testset "Volume" begin include("vol.jl") end
