using Test
using Crafts, Roly
using ForwardDiff, LinearAlgebra, NonlinearSolve, StaticArrays, Random, Statistics
using Convex, Clarabel
using JuMP, HiGHS
using Polymake

@testset "Crafts" verbose=true begin
    include("bondpotential.jl")
    include("assemblysystem.jl")
    include("yieldcalc.jl")
    include("reactions.jl")
    include("kinetics.jl")
    include("stability.jl")
    include("rates.jl")
    include("polyhedra.jl")
    include("facelattice.jl")
    include("designcone.jl")
    include("environmentrelations.jl")
    include("design.jl")
end