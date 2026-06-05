"""
    StructureCollection

A container type holding the most relevant physical attributes of collection of self-assembled structures, including the 
composition matrix `M`, the entropic partition functions `Zs`, and the aggregation `kernel`.
A `StructureCollection` can be constructed from an underlying `AssemblySystem`.
"""
struct StructureCollection{KF,P}
    M::Matrix{Int}
    Zs::Vector{Float64}
    kernel::KF
    sys::AssemblySystem
    strs::Vector{P}
end

Base.show(io::Core.IO, strs::StructureCollection) = print(io, "StructureCollection containing $(size(strs.M, 1)) structures")

"""
    StructureCollection(sys::AssemblySystem; veff=1e-2, enum_kwargs...)

Construct a `StructureCollection` from the binding rules and particle types defined by `sys`. 

This is done by enumerating all possible structures (use the optional `enum_kwargs` to modify the enumeration),
and then computing various attributes of the structure collection, employing **strong approximations**.

Computed attributes comprise:
- The composition matrix `M`
- The entropic partition functions `Zs`, which are approximated under the assumption of tree-like-ness, meaning that \$Z_s = \\|\\mathrm{SO}(d) / \\sigma_s v_\\mathrm{eff}\$.
- A trivial aggregation kernel that returns a rate of 1 for every reaction.
"""
function StructureCollection(sys::AssemblySystem; veff=1e-2, enum_kwargs...)
    strs = polygen(sys; enum_kwargs...)
    M = compositions(strs, sys)
    Zs = 8π^2 * inv.(s.σ for s in strs) / veff
    kernel(k, cut, components) = 1
    return StructureCollection(M, Zs, kernel, sys, strs)
end