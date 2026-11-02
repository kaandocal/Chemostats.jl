using Test
using Random
using Makie
using CairoMakie
using RecipesBase
using Chemostats

include("models/exponential.jl")

# Recipe defaults (`-->`) need a backend to declare supported keys; accept all here
RecipesBase.is_key_supported(::Symbol) = true

function simulated_chem(alg; L = 20, tmax = 3.0)
    Random.seed!(1)
    chem = Chemostat(make_population(L))
    Chemostats.simulate!(chem, tmax, alg; saveat = 0:0.5:tmax)
    chem
end

@testset "RecipesBase: Chemostat recipe" begin
    chem = simulated_chem(Chemostats.Strict(20))

    series = RecipesBase.apply_recipe(Dict{Symbol,Any}(), chem)
    subplots = [ rd.plotattributes[:subplot] for rd in series ]
    @test maximum(subplots) == 3

    xlabels = Dict(k => [ get(rd.plotattributes, :xlabel, nothing) for rd in series if rd.plotattributes[:subplot] == k ] for k in 1:3)
    @test all(isnothing, xlabels[1]) && all(isnothing, xlabels[2])
    @test xlabels[3][end] == "Time"

    ylabels = [ only(unique(rd.plotattributes[:ylabel] for rd in series if rd.plotattributes[:subplot] == k)) for k in 1:3 ]
    @test ylabels == ["log₁₀ N", "simulated size L", "growth rate Λ"]
end

@testset "RecipesBase: Lax adds thinning-rate panel" begin
    alg = Chemostats.Lax(20, 0.5; L_min = 5)
    chem = simulated_chem(alg)

    series = RecipesBase.apply_recipe(Dict{Symbol,Any}(:alg => alg), chem)
    @test maximum(rd.plotattributes[:subplot] for rd in series) == 4

    δ_series = [ rd for rd in series if rd.plotattributes[:subplot] == 4 ]
    @test only(unique(rd.plotattributes[:ylabel] for rd in δ_series)) == "thinning rate δ"
    @test get(δ_series[end].plotattributes, :xlabel, nothing) == "Time"
end

@testset "Makie: plot(chem)" begin
    fig = Makie.plot(simulated_chem(Chemostats.Strict(20)); Λ_gt = 1.0)
    axes = filter(b -> b isa Makie.Axis, fig.content)

    @test length(axes) == 3
    @test [ ax.xlabel[] for ax in axes ] == ["", "", "Time"]
    @test [ ax.ylabel[] for ax in axes ] == ["log₁₀ N", "simulated size L", "growth rate Λ"]
end

@testset "Makie: Lax adds thinning-rate panel" begin
    alg = Chemostats.Lax(20, 0.5; L_min = 5)
    fig = Makie.plot(simulated_chem(alg); Λ_gt = 1.0, alg = alg)
    axes = filter(b -> b isa Makie.Axis, fig.content)

    @test length(axes) == 4
    @test [ ax.xlabel[] for ax in axes ] == ["", "", "", "Time"]
    @test axes[end].ylabel[] == "thinning rate δ"
end
