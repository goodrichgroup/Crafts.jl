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
    return _relationsfrom(envs, rules, depth, nothing)
end

# Assemble the relations from an environment list. With `prescribed` bond classes the projection
# uses exactly those axes (in that order), so restricted relation sets stay composition-compatible
# with their parent.
function _relationsfrom(envs, rules, depth, prescribed)
    offset = Roly._markoffset(rules)
    ranges = _speciesranges(rules)
    ns = nspecies(rules)
    ne = length(envs)

    isempty(envs) && throw(ArgumentError("no environments to build relations from"))
    B = eltype(bondenvironments(first(envs)))
    rowidx = Dict{B,Int}()
    rowpairs = NTuple{2,B}[]
    entries = Dict{Tuple{Int,Int},Int}()
    clsidx = Dict{NTuple{2,Tuple{Int,Int}},Int}()
    bondclasses = NTuple{2,Tuple{Int,Int}}[]
    bondcounts = Dict{Tuple{Int,Int},Rational{Int}}()
    rootspecies = zeros(Int, ne)
    if prescribed !== nothing
        append!(bondclasses, prescribed)
        for (i, cls) in enumerate(prescribed)
            clsidx[cls] = i
        end
    end

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

    # order the bond axes deterministically; prescribed axes keep their given order
    if prescribed === nothing
        order = sortperm(bondclasses)
    else
        length(bondclasses) == length(prescribed) ||
            throw(ArgumentError("the environments carry a bond class the prescribed axes lack"))
        order = collect(eachindex(bondclasses))
    end
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

# The fiber of the composition direction `m`: `{μ ≥ 0 : relations⋅μ = 0, projection⋅μ = m}`.
# It is a polytope — every environment contributes +1 to its root-species row, so `Πμ = m`
# bounds `Σμ` — hence all LPs over it are bounded.
function _fibersystem(rel::EnvironmentRelations, m::AbstractVector)
    ne = length(rel.envs)
    nr = size(rel.relations, 1)
    d = size(rel.projection, 1)
    length(m) == d ||
        throw(DimensionMismatch("`m` has $(length(m)) entries but the composition space has $d"))
    A = [-Matrix{Rational{BigInt}}(I, ne, ne); Rational{BigInt}.(rel.relations);
         Rational{BigInt}.(rel.projection)]
    b = [zeros(Rational{BigInt}, ne + nr); Rational{BigInt}.(m)]
    return A, b, collect((ne + 1):(ne + nr + d))
end

"""
    realizablecounts(rel::EnvironmentRelations, m)

An exact count vector `μ ≥ 0` with `relations⋅μ = 0` and `projection⋅μ = m`, or `nothing` when
none exists. Feasibility certifies nothing about assembly; infeasibility rules `m` out exactly.
"""
function realizablecounts(rel::EnvironmentRelations, m::AbstractVector)
    A, b, linearity = _fibersystem(rel, m)
    status, x, _ = solvelp(A, b, zeros(Rational{Int}, length(rel.envs)); linearity)
    return status === :optimal ? x : nothing
end

"""
    fibersupport(rel::EnvironmentRelations, m)

Indices of every environment that can carry weight in a count vector realizing the composition
`m`: the support of the fiber `{μ ≥ 0 : relations⋅μ = 0, projection⋅μ = m}`. Empty when `m` is
infeasible.

Cropping maps any deeper fiber over `m` into this one, so only refinements of these environments
can matter when testing whether `m` survives at a larger depth — see [`refinementrelations`](@ref).
"""
function fibersupport(rel::EnvironmentRelations, m::AbstractVector)
    ne = length(rel.envs)
    A, b, linearity = _fibersystem(rel, m)
    w = zeros(Rational{Int}, ne)
    support = falses(ne)
    undecided = trues(ne)
    function absorb!(x)
        for j in 1:ne
            iszero(x[j]) && continue
            support[j] = true
            undecided[j] = false
        end
    end

    status, x, _ = solvelp(A, b, w; linearity)
    status === :optimal || return Int[]
    absorb!(x)
    for j in 1:ne
        undecided[j] || continue
        w[j] = 1
        status, x, _ = solvelp(A, b, w; linearity)
        w[j] = 0
        status === :optimal || error("the fiber LP must be bounded; solver returned `$status`")
        x[j] > 0 ? absorb!(x) : (undecided[j] = false)
    end
    return findall(support)
end

"""
    refinementrelations(rel::EnvironmentRelations, keep; depth=rel.depth+1, kwargs...)

The [`EnvironmentRelations`](@ref) of `rel.rules` at radius `depth`, restricted to environments
whose radius-`rel.depth` crop is one of `rel.envs[keep]`. The projection reuses `rel`'s bond-class
axes, so compositions are directly comparable across the two depths.

The restriction is lossless for a direction `m` whenever `keep` covers
[`fibersupport`](@ref)`(rel, m)`: any deeper count vector realizing `m` crops into the fiber, so
its support lies among these refinements. `realizablecounts(refined, m) === nothing` then rules
`m` out of the depth-`depth` outer cone exactly.
"""
function refinementrelations(rel::EnvironmentRelations, keep; depth=rel.depth + 1, kwargs...)
    kept = Set(rel.envs[collect(keep)])
    envs = eltype(rel.envs)[]
    particleenvironments((e, _) -> (crop(e, rel.depth) in kept && push!(envs, e); true),
                         rel.rules; depth, kwargs...)
    isempty(envs) && throw(ArgumentError("no refinements found for the kept environments"))
    return _relationsfrom(envs, rel.rules, depth, rel.bondclasses)
end

# Bond type index -> position of its symmetry class among `rel.bondclasses`.
function _bondclassmap(rel::EnvironmentRelations)
    rules = rel.rules
    sitelabel(loc) = labels(graphrep(Roly.species(rules, loc[1])))[first(
        bindingsites(Roly.species(rules, loc[1]), loc[2]).vertices)]
    return map(Roly.bonded_sites(rules)) do (locs1, locs2)
        l1, l2 = first(locs1), first(locs2)
        cls = minmax((l1[1], sitelabel(l1)), (l2[1], sitelabel(l2)))
        findfirst(==(cls), rel.bondclasses)
    end
end

# Composition of `poly` in `rel`'s symmetrized axes: species counts, then bond counts grouped by
# symmetry-equivalent bond type.
function _symcomposition(rel::EnvironmentRelations, classmap, poly)
    comp = composition(poly)
    ns = nspecies(rel.rules)
    c = zeros(Int, ns + length(rel.bondclasses))
    c[1:ns] .= comp[1:ns]
    for (β, cls) in enumerate(classmap)
        c[ns + cls] += comp[ns + β]
    end
    return c
end

"""
    findraywitness(rel::EnvironmentRelations, m; maxscale=1, maxsize=Inf)

Search for a finite structure whose symmetrized composition is proportional to `m`, by reverse
search pruned to compositions elementwise at most `maxscale` times the primitive representative
of `m`. Returns the witness `Polyform`, or `nothing` if the pruned search is exhausted.

A witness proves `m` realizable. Absence proves nothing: limit rays (infinite chains, lattices)
have no finite witness and need a periodic ansatz instead.
"""
function findraywitness(rel::EnvironmentRelations, m::AbstractVector; maxscale::Integer=1,
                        maxsize=Inf)
    d = size(rel.projection, 1)
    length(m) == d ||
        throw(DimensionMismatch("`m` has $(length(m)) entries but the composition space has $d"))
    mq = Rational{BigInt}.(m)
    mi = numerator.(mq .* lcm(denominator.(mq)))
    mi = mi .÷ reduce(gcd, mi)
    bound = maxscale .* mi

    classmap = _bondclassmap(rel)
    hit = Ref{Any}(nothing)
    function f(s, _)
        c = _symcomposition(rel, classmap, s)
        all(c .<= bound) || return REJECT
        if !iszero(c) && all(c[i] * mi[j] == c[j] * mi[i] for i in 1:d for j in (i + 1):d)
            hit[] = copy(s)
            return BREAK
        end
        return ACCEPT
    end
    polyenum(f, rel.rules; maxsize)
    return hit[]
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
