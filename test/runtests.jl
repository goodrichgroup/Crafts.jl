using Test
using Crafts, Roly
using ForwardDiff, LinearAlgebra, NonlinearSolve, StaticArrays, Random, Statistics
using Convex, Clarabel

@testset "Crafts" verbose=true begin
    include("bondpotential.jl")
    include("assemblysystem.jl")
    include("yieldcalc.jl")
    include("reactions.jl")
    include("kinetics.jl")
    include("stability.jl")
    include("rates.jl")
    include("lrs.jl")
    include("facelattice.jl")
    include("designcone.jl")
    include("design.jl")
end