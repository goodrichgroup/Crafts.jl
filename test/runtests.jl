using Test
using Crafts, Roly
using ForwardDiff, LinearAlgebra, NonlinearSolve, StaticArrays, Random, Statistics

@testset "Crafts" verbose=true begin
    include("bondpotential.jl")
    include("assemblysystem.jl")
    include("yieldcalc.jl")
    include("reactions.jl")
    include("kinetics.jl")
    include("stability.jl")
    include("rates.jl")
end