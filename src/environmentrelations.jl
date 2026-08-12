"""
    EnvironmentRelations

The exact linear relations satisfied by the particle-environment counts of every finite assembly,
at a fixed environment radius.

  - `envs`: every particle environment the rules allow; the count vector `μ` is indexed by it
  - `rowpairs`: the `(e, reverse(e))` bond-environment pair behind each relation
  - `relations`: `relations * μ == 0` exactly for the counts of any finite assembly or soup
  - `projection`: `projection * μ == [n_α; b̃]`, the composition with bond counts grouped by
    symmetry-equivalent bond type
  - `bondclasses`: the `((species, label), (species, label))` pair labelling each `b̃` axis
  - `depth`: environment radius

Each undirected bond is seen twice, once from each endpoint, and both sightings crop to the same
bond environment with the endpoint order swapped — so tail-counts of `e` and of `reverse(e)`
enumerate the same bonds, giving one exact integer relation per non-symmetric pair.
"""
struct EnvironmentRelations{E,B,S}
    envs::Vector{E}
    rowpairs::Vector{NTuple{2,B}}
    relations::SparseMatrixCSC{Int,Int}
    projection::Matrix{Rational{Int}}
    bondclasses::Vector{NTuple{2,Tuple{Int,Int}}}
    depth::Int
    rules::S
end

Base.show(io::IO, rel::EnvironmentRelations) =
    print(io, "EnvironmentRelations[k=", rel.depth, ", ", length(rel.envs), " environments, ",
          size(rel.relations, 1), " relations, d=", size(rel.projection, 1), "]")

_speciesranges(sys) = [extrema(labels(graphrep(Roly.species(sys, i)))) for i in 1:nspecies(sys)]
_speciesoflabel(ranges, l) = findfirst(r -> r[1] <= l <= r[2], ranges)

# The (species, label) description of a marked root vertex, stripping the marks.
function _rootclass(env, i, ranges, offset)
    l = mod1(labels(env.graph)[env.rootvertices[i]], offset)
    return (_speciesoflabel(ranges, l), l)
end

# Bond types up to particle symmetry: sites with equal labels are equivalent by Roly's label
# semantics, so the unordered pair of (species, label) is the finest bond resolution that is a
# class function of environments.
_bondclass(e, ranges, offset) = minmax(_rootclass(e, 1, ranges, offset),
                                       _rootclass(e, 2, ranges, offset))

"""
    environmentrelations(rules::BindingRules; depth=1, kwargs...)

Build the [`EnvironmentRelations`](@ref) of `rules` at environment radius `depth`. Additional
keyword arguments are passed to `particleenvironments`.
"""
function environmentrelations(rules::BindingRules; depth=1, kwargs...)
    envs = particleenvironments(rules; depth, kwargs...)
    offset = Roly._markoffset(rules)
    ranges = _speciesranges(rules)
    ns = nspecies(rules)
    ne = length(envs)

    B = eltype(bondenvironments(first(envs)))
    rowidx = Dict{B,Int}()
    rowpairs = NTuple{2,B}[]
    entries = Dict{Tuple{Int,Int},Int}()
    clsidx = Dict{NTuple{2,Tuple{Int,Int}},Int}()
    bondclasses = NTuple{2,Tuple{Int,Int}}[]
    bondcounts = Dict{Tuple{Int,Int},Rational{Int}}()
    rootspecies = zeros(Int, ne)

    for (j, w) in enumerate(envs)
        rootspecies[j] = _rootclass(w, 1, ranges, offset)[1]
        for e in bondenvironments(w)
            cls = _bondclass(e, ranges, offset)
            c = get!(clsidx, cls) do
                push!(bondclasses, cls)
                length(bondclasses)
            end
            bondcounts[(c, j)] = get(bondcounts, (c, j), 0//1) + 1//2

            ē = reverse(e)
            e == ē && continue   # a symmetric view balances itself; the relation would read 0 = 0
            if haskey(rowidx, e)
                r, s = rowidx[e], 1
            elseif haskey(rowidx, ē)
                r, s = rowidx[ē], -1
            else
                push!(rowpairs, (e, ē))
                rowidx[e] = length(rowpairs)
                r, s = length(rowpairs), 1
            end
            entries[(r, j)] = get(entries, (r, j), 0) + s
        end
    end

    # a particle has at most 4 bonds, so each column holds at most 4 nonzeros
    relations = sparse([r for (r, _) in keys(entries)], [j for (_, j) in keys(entries)],
                       collect(values(entries)), length(rowpairs), ne)
    dropzeros!(relations)

    # order the bond axes deterministically
    order = sortperm(bondclasses)
    projection = zeros(Rational{Int}, ns + length(bondclasses), ne)
    for j in 1:ne
        projection[rootspecies[j], j] = 1
    end
    for ((c, j), v) in bondcounts
        projection[ns + findfirst(==(c), order), j] = v
    end

    return EnvironmentRelations(envs, rowpairs, relations, projection, bondclasses[order],
                                Int(depth), rules)
end

"""
    environmentcounts(rel::EnvironmentRelations, poly::Polyform)

Count the particle environments of `poly` as a vector `μ` indexed like `rel.envs`.

Every environment of `poly` must be one the relations know; this holds whenever they were built
without truncation.
"""
function environmentcounts(rel::EnvironmentRelations, poly::Polyform)
    idx = Dict(e => j for (j, e) in enumerate(rel.envs))
    μ = zeros(Int, length(rel.envs))
    for e in particleenvironments(poly; depth=rel.depth)
        j = get(idx, e, 0)
        j == 0 && throw(ArgumentError("`poly` contains an environment the relations do not know"))
        μ[j] += 1
    end
    return μ
end

"""
    OuterCone

Outer bound on the composition cone: every composition of a finite assembly or soup satisfies
`facets * m ≤ 0`. Membership certifies nothing about realizability — the bound is one-sided.

  - `facets`: one inequality per row
  - `rays`: extreme rays, one per row

Entries are `Rational{Int}` whenever they fit, `Rational{BigInt}` otherwise (facet normals of
deep cones can carry huge coefficients).
"""
struct OuterCone{T<:Rational}
    facets::Matrix{T}
    rays::Matrix{T}
end

# Narrow to machine integers when possible; deep hulls can overflow them, exactness must not.
function _narrowcone(facets, rays)
    try
        return OuterCone(Rational{Int}.(facets), Rational{Int}.(rays))
    catch e
        e isa InexactError || rethrow()
        return OuterCone(Rational{BigInt}.(facets), Rational{BigInt}.(rays))
    end
end

Base.show(io::IO, c::OuterCone) =
    print(io, "OuterCone[d=", size(c.facets, 2), ", ", size(c.facets, 1), " facets, ",
          size(c.rays, 1), " rays]")

"""
    outercone(rel::EnvironmentRelations; method=:auto, optimizer=nothing)

Compute the outer composition cone `O = Π(𝒦)` exactly, where `𝒦 = {μ ≥ 0 : relations ⋅ μ = 0}`
is the feasible cone of the relations.

  - `method`: `:rays` enumerates the extreme rays of `𝒦` and projects them — exact but hopeless
    once `𝒦` has many extreme rays; `:project` gift-wraps `O` directly with one LP per face
    certificate, so its cost scales with the size of `O`, not of `𝒦`; `:auto` picks by size
  - `optimizer`: a JuMP-compatible optimizer (e.g. `HiGHS.Optimizer`; requires JuMP loaded) to
    solve the `:project` LPs in floating point with rationalized results — much faster, but the
    cone is only as exact as the rationalization

Soundness: the counts of every finite assembly lie in `𝒦`, so every composition lies in `O`.
Feasible points of `𝒦` are not claimed realizable; only the outer bound is exact.
"""
function outercone(rel::EnvironmentRelations; method::Symbol=:auto, optimizer=nothing)
    ne = length(rel.envs)
    nr = size(rel.relations, 1)
    if method === :auto
        method = optimizer === nothing && ne <= 64 ? :rays : :project
    end
    linearity = collect((ne + 1):(ne + nr))

    if method === :rays
        A = vcat(-Matrix{Rational{Int}}(I, ne, ne), Rational{Int}.(rel.relations))
        rays, _ = extremerays(A, zeros(Rational{Int}, ne + nr); linearity)
        projected = Rational{Int}.(rays) * transpose(rel.projection)
        facets, bf, _ = facetsof(projected; rays=true)
        minimalrays, _ = extremerays(Rational{Int}.(facets), Rational{Int}.(bf))
        return _narrowcone(facets, minimalrays)
    elseif method === :project
        # environments whose relation column vanishes are feasible on their own (the monomer
        # environments always are); their projections seed the hull for free
        # sparse: at large nₑ the dense identity block alone would be unallocatable
        A = [-sparse(one(Rational{Int}) * I, ne, ne); sparse(Rational{Int}.(rel.relations))]
        seedcols = findall(j -> nnz(rel.relations[:, j]) == 0, 1:ne)
        seeds = permutedims(rel.projection[:, seedcols])
        facets, rays = projectcone(rel.projection, A; linearity, seeds, optimizer)
        return _narrowcone(facets, rays)
    end
    throw(ArgumentError("unknown method `$(repr(method))`; expected `:auto`, `:rays` or `:project`"))
end
