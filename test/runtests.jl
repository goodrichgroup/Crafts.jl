using Test
using Crafts, Roly
using ForwardDiff, LinearAlgebra, NonlinearSolve

@testset "Crafts" verbose=true begin
    include("assemblysystem.jl")
    include("yieldcalc.jl")
    include("reactions.jl")
end