"""
    FaceLattice

The faces of a cone or polytope, grouped by rank.

  - `levels[k+1]`: the faces of rank `k`, each a `BitVector` over the generators; rank 0 is the
    smallest face, the apex of a pointed cone or the empty face of a polytope
  - `dualsets[k+1]`: the facets containing each face, when `facetsets=true` was passed
  - `covers[k+1]`: pairs `(i, j)` linking face `i` of rank `k` to face `j` of rank `k+1`,
    when `covers=true` was passed
  - `dual`: whether the lattice was built from the facets rather than the generators
"""
struct FaceLattice
    levels::Vector{Vector{BitVector}}
    dualsets::Union{Nothing,Vector{Vector{BitVector}}}
    covers::Union{Nothing,Vector{Vector{Tuple{Int,Int}}}}
    dual::Bool
end

"""
    faces(L::FaceLattice, k)

Return the faces of rank `k`, each a `BitVector` over the generators.
"""
faces(L::FaceLattice, k::Integer) = L.levels[k + 1]

"""
    fvector(L::FaceLattice)

Return the number of faces of each rank, starting at rank zero.
"""
fvector(L::FaceLattice) = length.(L.levels)

"""
    covers(L::FaceLattice, k)

Return the cover relations from rank `k` to `k+1` as `(i, j)` index pairs.
"""
function covers(L::FaceLattice, k::Integer)
    L.covers === nothing &&
        throw(ArgumentError("this lattice was built without cover relations; pass `covers=true`"))
    return L.covers[k + 1]
end

Base.show(io::IO, L::FaceLattice) =
    print(io, "FaceLattice[f=", join(fvector(L), ","), L.dual ? ", dual]" : "]")

"""
    facelattice(incidence; maxdim, dual=false, covers=false, facetsets=false)

Enumerate the faces of a cone or polytope from its generator-facet incidences.

  - `incidence[i, j]`: whether generator `j` lies on facet `i`
  - `maxdim`: stop after this rank; the default enumerates the whole lattice
  - `dual`: walk down from the facets instead of up from the generators
  - `covers`: also return the cover relations
  - `facetsets`: also tag each face with the facets containing it

Returns a [`FaceLattice`](@ref). Requires Polymake to be loaded.
"""
function facelattice(incidence::AbstractMatrix{Bool}; maxdim=typemax(Int), dual=false,
                     covers=false, facetsets=false)
    inc = dual ? permutedims(incidence) : incidence
    levels, duals, edges = _facelattice(_backend(), Matrix{Bool}(inc), maxdim, covers, facetsets)
    return FaceLattice(levels, duals, edges, dual)
end

"""
    facelattice(A, b; kwargs...)

Enumerate the faces of the cone `A*x <= b`, computing its generators and incidences first.
"""
function facelattice(A::AbstractMatrix{<:ExactNum}, b::AbstractVector{<:ExactNum}; linearity=Int[],
                     kwargs...)
    rays, verts = extremerays(A, b; linearity)
    V = vcat(rays, verts)
    isempty(V) && throw(ArgumentError("the system `A*x <= b` has no generators"))
    flags = [trues(size(rays, 1)); falses(size(verts, 1))]
    _, _, inc = facetsof(V; rays=flags)
    return facelattice(inc; kwargs...)
end
