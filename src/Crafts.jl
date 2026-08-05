module Crafts

using LinearAlgebra, ForwardDiff, Graphs, Statistics, NonlinearSolve, LogExpFunctions, OrdinaryDiffEq, StaticArrays
using Roly, NautyGraphs
import Roly: nspecies, nbonds, dimension


export logdensities, densities, logyields, yields, logparticledensities, particledensities, chemicalpotentials
export density_jacobian, yield_jacobian, particledensity_jacobian

export make_harmonicpotential, map_potential, hessian, entropy, EntropySolver, EntropyModel
export TetheredLaplace, COMLaplace, TreeApproximation
export AssemblySystem, structures, nstructures, countstructures, iscomplete, compositionmatrix, partitionfunctions, simulate_kinetics

include("utils.jl")
include("bondpotentials.jl")
include("entropy.jl")
include("assemblysystem.jl")
include("yieldcalc.jl")
include("kinetics.jl")

end # module Crafts
