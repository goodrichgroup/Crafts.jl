"""
    parametermap(asys::AssemblySystem; bondgroups=nothing)

Build the matrix `P` relating a reduced set of parameters `η` to the full set, `ξ = P*η`.

Experiments often cannot tune every bond energy independently. Grouping bond types that share an
energy shrinks the design space, which is usually what makes the constraint cone small enough to
analyse.

  - `bondgroups`: a partition of the bond types into groups that share an energy, or `:uniform` to
    put them all in one group. `nothing` leaves every bond independent
"""
function parametermap(asys::AssemblySystem; bondgroups=nothing)
    nμ, nb = nspecies(asys), nbonds(asys)
    bondgroups === nothing && return Matrix{Int}(I, nμ + nb, nμ + nb)

    groups = bondgroups === :uniform ? [collect(1:nb)] : [collect(Int, g) for g in bondgroups]
    flat = reduce(vcat, groups; init=Int[])
    sort(flat) == 1:nb ||
        throw(ArgumentError("`bondgroups` must partition the $nb bond types, got $(sort(flat))"))

    P = zeros(Int, nμ + nb, nμ + length(groups))
    for i in 1:nμ
        P[i, i] = 1
    end
    for (g, bs) in pairs(groups), bond in bs
        P[nμ + bond, nμ + g] = 1
    end
    return P
end

_reducedM(asys; bondgroups=nothing, kwargs...) =
    compositionmatrix(asys; kwargs...) * parametermap(asys; bondgroups)

"""
    ConstraintCone

The thermodynamic constraint cone `{ξ : Mₛ·ξ ≤ 0 for every structure s}` of an assembly system.

Directions in the cone are the ways the parameters can be pushed to their asymptotic limit. A
direction lying on the boundary assembles exactly those structures whose constraint is tight there,
so the faces of the cone enumerate the sets of structures that can be assembled together.

  - `M`: the composition matrix, one constraint per structure
  - `designable`: structures whose constraint is not redundant, i.e. designable on their own
  - `rays`: the extreme rays of the cone, one per row
  - `incidence[s, r]`: whether structure `s` is tight on ray `r`, for every structure

Incidence covers all structures, not just the designable ones. A redundant constraint still touches
the cone on some lower-dimensional face, and those touchings are what make a structure assemble
alongside others rather than alone.
"""
struct ConstraintCone
    M::Matrix{Rational{BigInt}}
    designable::Vector{Int}
    rays::Matrix{Rational{BigInt}}
    incidence::BitMatrix
end

Base.show(io::IO, c::ConstraintCone) =
    print(io, "ConstraintCone[d=", size(c.M, 2), ", ", size(c.M, 1), " structures, ",
          length(c.designable), " designable, ", size(c.rays, 1), " rays]")

"""
    constraintcone(asys::AssemblySystem; kwargs...)

Build the thermodynamic constraint cone of `asys`.

This enumerates the extreme rays, which can be expensive; [`designablestructures`](@ref) is much
cheaper when only the individually designable structures are needed.
"""
function constraintcone(asys::AssemblySystem; bondgroups=nothing, kwargs...)
    M = Rational{BigInt}.(_reducedM(asys; bondgroups, kwargs...))
    b = zeros(Rational{BigInt}, size(M, 1))

    _, _, designable = removeredundancy(M, b)
    rays, _ = extremerays(M, b)

    # Incidence is taken against the original rows rather than round-tripping through `facetsof`,
    # which would reorder the facets and lose track of which structure each one came from.
    inc = falses(size(M, 1), size(rays, 1))
    for s in axes(M, 1), r in axes(rays, 1)
        inc[s, r] = iszero(sum(M[s, k] * rays[r, k] for k in axes(rays, 2)))
    end

    # The cone is always pointed: the monomers supply every species direction, and subtracting them
    # from a dimer supplies every bond direction, so `M` has full rank and there is no lineality.
    return ConstraintCone(M, designable, rays, inc)
end

"""
    designablestructures(asys::AssemblySystem; kwargs...)

Return the structures that can each be assembled at 100% yield on their own.

These are the structures whose thermodynamic constraint is not redundant. Finding them needs only
redundancy removal, so this stays affordable on systems where the full cone does not.
"""
function designablestructures(asys::AssemblySystem; bondgroups=nothing, kwargs...)
    M = Rational{BigInt}.(_reducedM(asys; bondgroups, kwargs...))
    _, _, kept = removeredundancy(M, zeros(Rational{BigInt}, size(M, 1)))
    return kept
end

designablestructures(c::ConstraintCone) = c.designable

_asidxs(idxs) = idxs isa Number ? [Int(idxs)] : sort!(collect(Int, idxs))

"""
    minimaldesignableset(c::ConstraintCone, idxs)

Return the smallest set of structures that is designable and contains `idxs`.

Designable sets are closed under intersection, so this set is unique. Assembling a structure that is
not designable on its own means assembling its minimal designable set, and then trading off the
yields within it.
"""
function minimaldesignableset(c::ConstraintCone, idxs)
    want = _asidxs(idxs)
    all(in(axes(c.M, 1)), want) || throw(ArgumentError("structure index out of range"))

    # The face on which every structure in `want` is tight is spanned by the rays where they all
    # vanish. Closing that set means adding every structure tight on all of those rays.
    onface = trues(size(c.rays, 1))
    for s in want
        onface .&= view(c.incidence, s, :)
    end
    # No ray survives, so the only such point is the apex, where every structure has density Ωₛ.
    any(onface) || return collect(axes(c.M, 1))

    return _activeset(c, onface)
end

# The structures tight on every ray of a face, which is the designable set that face carries.
function _activeset(c::ConstraintCone, onface::AbstractVector{Bool})
    active = Int[]
    for s in axes(c.M, 1)
        all(r -> !onface[r] || c.incidence[s, r], axes(c.rays, 1)) && push!(active, s)
    end
    return active
end

"""
    isdesignable(c::ConstraintCone, idxs)

Whether the structures `idxs` can be assembled together at a combined 100% yield.
"""
isdesignable(c::ConstraintCone, idxs) = minimaldesignableset(c, idxs) == _asidxs(idxs)

"""
    necessarychimeras(c::ConstraintCone, idxs)

Return the structures that cannot be suppressed while assembling `idxs`.

These are the off-target structures that any parameter choice must tolerate, so they set the ceiling
on the achievable yield of the targets.
"""
necessarychimeras(c::ConstraintCone, idxs) = setdiff(minimaldesignableset(c, idxs), _asidxs(idxs))

"""
    designablesets(c::ConstraintCone; maxsize=typemax(Int))

Return every designable set of at most `maxsize` structures, smallest first.

Small sets sit on the high-dimensional faces of the cone, which are the cheap end of the face
lattice, so `maxsize` keeps this affordable on systems whose full lattice is out of reach. The
cone's interior, where nothing assembles at all, is not listed.
"""
function designablesets(c::ConstraintCone; maxsize=typemax(Int))
    # Only the non-redundant constraints bound the cone, so the lattice is built from those. Walking
    # it from the facets down means the small sets, which are the interesting ones, come first.
    L = facelattice(c.incidence[c.designable, :]; dual=true, maxdim=maxsize, facetsets=true)

    # A face is identified by the rays it spans; its designable set is every structure tight on all
    # of them, redundant constraints included.
    sets = Vector{Int}[]
    for level in L.dualsets, onface in level
        active = _activeset(c, onface)
        isempty(active) || length(active) > maxsize || push!(sets, active)
    end
    unique!(sets)
    return sort!(sets; by=length)
end

"""
    codimension(asys::AssemblySystem, idxs)

Return the number of parameter directions that change the densities within the set `idxs`.

This is the codimension `c_f` of the face the set sits on. It counts the dimensions of the reduced
design space left once the parameters are aligned with that face.
"""
function codimension(asys::AssemblySystem, idxs; bondgroups=nothing)
    return rank(Float64.(_reducedM(asys; bondgroups)[_asidxs(idxs), :]))
end

codimension(c::ConstraintCone, idxs) = rank(Float64.(c.M[_asidxs(idxs), :]))

# Whether some parameter direction raises every density in the set by the same factor, which changes
# the total amount assembled but not the yields relative to each other.
function _hasuniformdirection(MT)
    n = size(MT, 1)
    return rank(hcat(MT, ones(n))) == rank(MT)
end

"""
    relativeyielddofs(asys::AssemblySystem, idxs)

Return how many independent parameters control the yields of `idxs` relative to each other.

This is [`codimension`](@ref) less one whenever some direction rescales every density in the set
equally, since that leaves the relative yields untouched.
"""
function relativeyielddofs(asys::AssemblySystem, idxs; bondgroups=nothing)
    MT = Float64.(_reducedM(asys; bondgroups)[_asidxs(idxs), :])
    return rank(MT) - (_hasuniformdirection(MT) ? 1 : 0)
end

"""
    yielddirections(asys::AssemblySystem, idxs)

Return an orthonormal basis for the parameter directions that tune the yields within `idxs`.

Columns span the reduced design space of the set, with any direction that merely rescales all of its
densities removed. The number of columns is [`relativeyielddofs`](@ref).
"""
function yielddirections(asys::AssemblySystem, idxs; bondgroups=nothing)
    MT = Float64.(_reducedM(asys; bondgroups)[_asidxs(idxs), :])
    r = rank(MT)
    r == 0 && return zeros(size(MT, 2), 0)

    F = svd(MT)
    basis = F.V[:, 1:r]                     # the row space: everything else leaves the set untouched

    if _hasuniformdirection(MT)
        u = pinv(MT) * ones(size(MT, 1))    # raises every density equally
        n = norm(u)
        if n > 0
            u ./= n
            basis = basis - u * (u' * basis)
            basis = Matrix(qr(basis).Q)[:, 1:(r - 1)]
        end
    end
    return basis
end

"""
    relativeyields(asys::AssemblySystem, idxs, ζ)

Return the yields of `idxs` relative to each other at reduced-design-space coordinates `ζ`.

Sweeping `ζ` traces out every combination of relative yields the set can reach, whatever the binding
energies and concentrations. `ζ` has one entry per column of [`yielddirections`](@ref).

  - returns the yields of `idxs` normalized to sum to one
"""
function relativeyields(asys::AssemblySystem, idxs, ζ; bondgroups=nothing)
    want = _asidxs(idxs)
    Y = yielddirections(asys, want; bondgroups)
    length(ζ) == size(Y, 2) ||
        throw(DimensionMismatch("`ζ` has $(length(ζ)) entries but the set has $(size(Y, 2)) " *
                                "relative-yield directions"))

    logρ = log.(partitionfunctions(asys)[want]) .+
           _reducedM(asys; bondgroups)[want, :] * (Y * collect(ζ))
    return exp.(logρ .- logsumexp(logρ))
end
