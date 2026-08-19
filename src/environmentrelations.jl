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
    slacks::SparseMatrixCSC{Int,Int}      # inequality rows `slacks * μ ≤ 0` (wildcard tier)
    projection::Matrix{Rational{Int}}
    bondclasses::Vector{NTuple{2,Tuple{Int,Int}}}
    depth::Int
    rules::S
end

Base.show(io::IO, rel::EnvironmentRelations) =
    print(io, "EnvironmentRelations[k=", rel.depth, ", ", length(rel.envs), " environments, ",
          size(rel.relations, 1), " relations",
          size(rel.slacks, 1) == 0 ? "" : ", $(size(rel.slacks, 1)) slack rows",
          ", d=", size(rel.projection, 1), "]")

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
# with their parent. `resolve(w, e)` maps each bond crop to the resolution its balance row lives
# at — the identity for uniform depth, the coarsest-common rule for adaptive depth.
function _relationsfrom(envs, rules, depth, prescribed, resolve=(w, e) -> e)
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
        for e0 in bondenvironments(w)
            cls = _bondclass(e0, ranges, offset)
            c = get!(clsidx, cls) do
                push!(bondclasses, cls)
                length(bondclasses)
            end
            bondcounts[(c, j)] = get(bondcounts, (c, j), 0//1) + 1//2

            e = resolve(w, e0)
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

    return EnvironmentRelations(envs, rowpairs, relations, spzeros(Int, 0, ne), projection,
                                bondclasses[order], Int(depth), rules)
end

"""
    environmentcounts(rel::EnvironmentRelations, poly::Polyform)

Count the particles of `poly` into `rel`'s environment bins, as a vector `μ` indexed like
`rel.envs`. Each particle lands in its deepest matching bin, so uniform and adaptive-depth
relation systems are both covered.

Every particle must match some bin; this holds whenever the relations were built without
truncation.
"""
function environmentcounts(rel::EnvironmentRelations, poly::Polyform)
    idx = Dict(e => j for (j, e) in enumerate(rel.envs))
    depths = sort!(unique(e.depth for e in rel.envs); rev=true)
    μ = zeros(Int, length(rel.envs))
    for p in 1:nparticles(poly)
        j = 0
        for k in depths
            j = get(idx, ParticleEnvironment(poly, p; depth=k), 0)
            j == 0 || break
        end
        j == 0 &&
            throw(ArgumentError("`poly` contains an environment the relations do not know"))
        μ[j] += 1
    end
    return μ
end

"""
    adaptiverelations(rel::EnvironmentRelations, refine; kwargs...)

The exact adaptive-depth relations: the environments `rel.envs[refine]` are replaced by their
complete refinement families one radius deeper, the rest stay as they are, and no slack is
needed anywhere.

Each bond's balance row lives at the coarsest common resolution of its two sides: a
refined-side environment reads its neighbor's coarse environment out of the bond crop
([`rootenvironment`](@extref)) and drops to the coarse pairing when the neighbor is unrefined.
The environment count interpolates between `rel`'s and the uniformly deeper system's, one
refinement family at a time. Additional keyword arguments are passed to `particleenvironments`.
"""
function adaptiverelations(rel::EnvironmentRelations, refine; kwargs...)
    _requireuniform(rel, "adaptiverelations")
    isempty(refine) && return rel
    refined = Set(rel.envs[collect(refine)])
    k = rel.depth

    bins = eltype(rel.envs)[e for e in rel.envs if !(e in refined)]
    nshallow = length(bins)
    particleenvironments((e, _) -> (crop(e, k) in refined && push!(bins, e); true),
                         rel.rules; depth=k + 1, kwargs...)
    length(bins) > nshallow ||
        throw(ArgumentError("the environments to refine have no refinements"))

    # coarsest-common pairing: the unrefined side's crops are already at the coarse resolution;
    # the refined side downgrades exactly when the neighbor read from the bond crop is unrefined
    resolve(w, e) = w.depth == k ? e :
                    rootenvironment(e, 2, k) in refined ? e : crop(e, k - 1)
    return _relationsfrom(bins, rel.rules, k + 1, rel.bondclasses, resolve)
end

# Standard balance-row bookkeeping: orient the (e, reverse(e)) pair, create the row on first
# sighting, accumulate ±1.
function _pairentry!(rowidx, rowpairs, entries, e, j)
    ē = reverse(e)
    e == ē && return
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
    return
end

"""
    wildcardrelations(rel::EnvironmentRelations, listed; kwargs...)

The wildcard-tier system for an arbitrary partial listing of one-radius-deeper refinements:
listed refinements become bins of their own; each coarse environment with unlisted refinements
stays as its own leftover bin, catching every particle no listed refinement matches.

Rows come in two kinds. Coarse balance rows stay exact equalities (every bin's coarse bond
crops are determined). Fine balance rows carry the listed bins' exact counts, with the leftover
bins' unknown fine sightings bounded by their coarse bond counts — a pair of valid inequalities
per fine bond-environment pair (`slacks`). With a complete listing the leftover bins vanish and
the system is equivalent to the uniformly deeper one; with an empty listing it degenerates to
`rel`. Sound for the full chemistry at any listing, tightening as the listing grows.
"""
function wildcardrelations(rel::EnvironmentRelations, listed; kwargs...)
    _requireuniform(rel, "wildcardrelations")
    k = rel.depth
    E = eltype(rel.envs)
    listedset = Set{E}(listed)
    all(e -> e.depth == k + 1, listedset) ||
        throw(ArgumentError("`listed` must hold depth-$(k + 1) environments"))

    bins = E[]
    familysize = Dict{E,Int}()
    listedcount = Dict{E,Int}()
    particleenvironments((e, _) -> begin
                             w = crop(e, k)
                             familysize[w] = get(familysize, w, 0) + 1
                             if e in listedset
                                 push!(bins, e)
                                 listedcount[w] = get(listedcount, w, 0) + 1
                             end
                             true
                         end, rel.rules; depth=k + 1, kwargs...)
    leftovers = E[w for w in rel.envs if get(listedcount, w, 0) < familysize[w]]
    envs = vcat(bins, leftovers)
    ne = length(envs)

    offset = Roly._markoffset(rel.rules)
    ranges = _speciesranges(rel.rules)
    ns = nspecies(rel.rules)
    B = eltype(bondenvironments(first(rel.envs)))

    crowidx = Dict{B,Int}()   # coarse equality rows
    crowpairs = NTuple{2,B}[]
    centries = Dict{Tuple{Int,Int},Int}()
    frowidx = Dict{B,Int}()   # fine rows, to be emitted with slack
    frowpairs = NTuple{2,B}[]
    fentries = Dict{Tuple{Int,Int},Int}()
    leftcoarse = Dict{Tuple{Int,B},Int}()   # leftover bin × coarse bond env -> bond count
    clsidx = Dict{NTuple{2,Tuple{Int,Int}},Int}()
    bondclasses = copy(rel.bondclasses)
    for (i, cls) in enumerate(bondclasses)
        clsidx[cls] = i
    end
    bondcounts = Dict{Tuple{Int,Int},Rational{Int}}()
    rootspecies = zeros(Int, ne)

    for (j, w) in enumerate(envs)
        rootspecies[j] = _rootclass(w, 1, ranges, offset)[1]
        fine = w.depth == k + 1
        for be in bondenvironments(w)
            cls = _bondclass(be, ranges, offset)
            haskey(clsidx, cls) ||
                throw(ArgumentError("the refinements carry a bond class the parent lacks"))
            c = clsidx[cls]
            bondcounts[(c, j)] = get(bondcounts, (c, j), 0//1) + 1//2

            ce = fine ? crop(be, k - 1) : be
            _pairentry!(crowidx, crowpairs, centries, ce, j)
            if fine
                _pairentry!(frowidx, frowpairs, fentries, be, j)
            else
                leftcoarse[(j, ce)] = get(leftcoarse, (j, ce), 0) + 1
            end
        end
    end

    relations = sparse([r for (r, _) in keys(centries)], [j for (_, j) in keys(centries)],
                       collect(values(centries)), length(crowpairs), ne)
    dropzeros!(relations)

    # two slack rows per fine pair: the listed counts' imbalance in either direction is bounded
    # by the leftover bins' matching coarse bond counts
    Is, Js, Vs = Int[], Int[], Int[]
    for ((r, j), v) in fentries
        push!(Is, 2r - 1); push!(Js, j); push!(Vs, v)
        push!(Is, 2r); push!(Js, j); push!(Vs, -v)
    end
    for (r, (e, ē)) in enumerate(frowpairs)
        ce, cē = crop(e, k - 1), crop(ē, k - 1)
        for (j, w) in enumerate(envs)
            w.depth == k || continue
            # Σ_listed (N_e − N_ē) μ = hidden_ē − hidden_e, so the positive side is bounded by
            # the leftovers' coarse sightings of ē, the negative side by those of e
            n⁺ = get(leftcoarse, (j, cē), 0)
            n⁺ == 0 || (push!(Is, 2r - 1); push!(Js, j); push!(Vs, -n⁺))
            n⁻ = get(leftcoarse, (j, ce), 0)
            n⁻ == 0 || (push!(Is, 2r); push!(Js, j); push!(Vs, -n⁻))
        end
    end
    slacks = sparse(Is, Js, Vs, 2 * length(frowpairs), ne)
    dropzeros!(slacks)

    projection = zeros(Rational{Int}, ns + length(bondclasses), ne)
    for j in 1:ne
        projection[rootspecies[j], j] = 1
    end
    for ((c, j), v) in bondcounts
        projection[ns + c, j] = v
    end

    return EnvironmentRelations(envs, crowpairs, relations, slacks, projection, bondclasses,
                                k + 1, rel.rules)
end

"""
    listingrelations(rel::EnvironmentRelations, budget; priority=Int[], closure=true, kwargs...)

One-knob coarseness control: refine within a bin budget. Refinement families are swapped in
whole while they fit (the exact adaptive tier); the first family that no longer fits is refined
partially up to the budget, with its remainder caught by a leftover bin (the wildcard tier).
`budget` counts total bins; `budget = length(rel.envs)` changes nothing, a large budget recovers
the uniformly deeper system.

  - `priority`: indices into `rel.envs` refined first — typically the union of the eliminated
    rays' fiber supports (`certifyrays` returns them per verdict). Refinement only bites where
    the damage is; an unguided budget can be spent entirely on irrelevant families
  - `closure`: also prioritize (second) the environments adjacent to a priority family across a
    bond, read from the refinements' bond crops — the fine balance rows that cut a fake ray need
    both sides of their bonds refined

Within each group (priority, adjacency closure, rest) smaller families go first.
"""
function listingrelations(rel::EnvironmentRelations, budget::Integer; priority=Int[],
                          closure::Bool=true, kwargs...)
    _requireuniform(rel, "listingrelations")
    budget >= length(rel.envs) ||
        throw(ArgumentError("`budget` is below the current bin count $(length(rel.envs))"))
    k = rel.depth
    E = eltype(rel.envs)
    fams = Dict{E,Vector{E}}()
    prioset = Set{E}(rel.envs[collect(priority)])
    nbrset = Set{E}()
    particleenvironments((e, _) -> begin
                             w = crop(e, k)
                             push!(get!(fams, w, E[]), e)
                             if closure && w in prioset
                                 for be in bondenvironments(e)
                                     push!(nbrset, rootenvironment(be, 2, k))
                                 end
                             end
                             true
                         end, rel.rules; depth=k + 1, kwargs...)

    group(w) = w in prioset ? 0 : w in nbrset ? 1 : 2
    order = sortperm([(group(w), length(fams[w])) for w in rel.envs])
    listed = E[]
    total = length(rel.envs)
    for i in order
        f = fams[rel.envs[i]]
        if total + length(f) - 1 <= budget
            append!(listed, f)
            total += length(f) - 1
        else
            m = budget - total
            m > 0 && append!(listed, f[1:m])
            break
        end
    end
    return wildcardrelations(rel, listed)
end
# It is a polytope — every environment contributes +1 to its root-species row, so `Πμ = m`
# bounds `Σμ` — hence all LPs over it are bounded.
function _fibersystem(rel::EnvironmentRelations, m::AbstractVector)
    ne = length(rel.envs)
    ns = size(rel.slacks, 1)
    nr = size(rel.relations, 1)
    d = size(rel.projection, 1)
    length(m) == d ||
        throw(DimensionMismatch("`m` has $(length(m)) entries but the composition space has $d"))
    A = [-Matrix{Rational{BigInt}}(I, ne, ne); Rational{BigInt}.(rel.slacks);
         Rational{BigInt}.(rel.relations); Rational{BigInt}.(rel.projection)]
    b = [zeros(Rational{BigInt}, ne + ns + nr); Rational{BigInt}.(m)]
    return A, b, collect((ne + ns + 1):(ne + ns + nr + d))
end

"""
    realizablecounts(rel::EnvironmentRelations, m; optimizer=nothing)

A count vector `μ ≥ 0` with `relations⋅μ = 0` and `projection⋅μ = m`, or `nothing` when none
exists. Feasibility certifies nothing about assembly; `nothing` rules `m` out exactly.

Without an `optimizer` the fiber LP is solved exactly through lrs — viable only for small
systems. With one (requires JuMP loaded), the LP is solved in floating point and an
infeasibility verdict is only reported after its Farkas certificate verifies in exact
arithmetic; an unverifiable certificate returns `missing` (unknown) instead, so `nothing`
always remains an exact elimination.
"""
function realizablecounts(rel::EnvironmentRelations, m::AbstractVector; optimizer=nothing)
    if optimizer !== nothing
        d = size(rel.projection, 1)
        length(m) == d ||
            throw(DimensionMismatch("`m` has $(length(m)) entries but the composition space has $d"))
        return _fiberfeasible(optimizer, rel.slacks, rel.relations, rel.projection, m)
    end
    length(rel.envs) <= 4096 ||
        throw(ArgumentError("the exact fiber LP is infeasible at $(length(rel.envs)) " *
                            "environments; pass a JuMP `optimizer`"))
    A, b, linearity = _fibersystem(rel, m)
    status, x, _ = solvelp(A, b, zeros(Rational{Int}, length(rel.envs)); linearity)
    return status === :optimal ? x : nothing
end

_fiberfeasible(optimizer, S, L, P, m) =
    throw(ArgumentError("passing an `optimizer` requires JuMP to be loaded, e.g. `using JuMP, HiGHS`"))

_requireuniform(rel::EnvironmentRelations, what) =
    (size(rel.slacks, 1) == 0 && allequal(e.depth for e in rel.envs)) ||
    throw(ArgumentError("$what requires a uniform, slack-free relation system"))

_wedgemembers(optimizer, L, P, h) =
    throw(ArgumentError("passing an `optimizer` requires JuMP to be loaded, e.g. `using JuMP, HiGHS`"))

_fibermembers(optimizer, S, L, P, m) =
    throw(ArgumentError("passing an `optimizer` requires JuMP to be loaded, e.g. `using JuMP, HiGHS`"))

"""
    fibersupport(rel::EnvironmentRelations, m; optimizer=nothing)

Indices of every environment that can carry weight in a count vector realizing the composition
`m`: the support of the fiber `{μ ≥ 0 : relations⋅μ = 0, projection⋅μ = m}`. Empty when `m` is
infeasible.

Cropping maps any deeper fiber over `m` into this one, so only refinements of these environments
can matter when testing whether `m` survives at a larger depth — see [`refinementrelations`](@ref).

Without an `optimizer` the scan runs one exact lrs LP per environment, which is one subprocess
each and dominates certification on any system of size. With one (requires JuMP loaded) the same
scan runs in floating point on a single model with the objective swapped per environment.
Over-collecting is sound here — the support only decides which environments get refined, so a
superset costs refinement work and never a verdict — hence float membership needs no certificate.
"""
function fibersupport(rel::EnvironmentRelations, m::AbstractVector; optimizer=nothing)
    optimizer === nothing ||
        return _fibermembers(optimizer, rel.slacks, rel.relations, rel.projection, m)
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

# Support of the wedge slab `{μ ≥ 0 : Lμ = 0, (hᵀΠ)μ = 1}`: every environment that can carry
# weight in a count vector whose composition violates the survivor-cone facet `h`. Unlike a
# fiber, the slab can be unbounded; an unbounded coordinate is a member.
function _wedgesupport(rel::EnvironmentRelations, h::AbstractVector)
    ne = length(rel.envs)
    ns = size(rel.slacks, 1)
    nr = size(rel.relations, 1)
    w = vec(transpose(Rational{BigInt}.(h)) * rel.projection)
    A = [-Matrix{Rational{BigInt}}(I, ne, ne); Rational{BigInt}.(rel.slacks);
         Rational{BigInt}.(rel.relations); permutedims(w)]
    b = [zeros(Rational{BigInt}, ne + ns + nr); one(Rational{BigInt})]
    linearity = collect((ne + ns + 1):(ne + ns + nr + 1))

    support = falses(ne)
    undecided = trues(ne)
    obj = zeros(Rational{Int}, ne)
    function absorb!(x)
        for j in 1:ne
            iszero(x[j]) && continue
            support[j] = true
            undecided[j] = false
        end
    end

    status, x, _ = solvelp(A, b, obj; linearity)
    status === :optimal || return Int[]   # empty slab: no feasible counts beyond this facet
    absorb!(x)
    for j in 1:ne
        undecided[j] || continue
        obj[j] = 1
        status, x, _ = solvelp(A, b, obj; linearity)
        obj[j] = 0
        if status === :unbounded
            support[j] = true
            undecided[j] = false
        elseif status === :optimal
            x[j] > 0 ? absorb!(x) : (undecided[j] = false)
        else
            error("the wedge LP must be feasible; solver returned `$status`")
        end
    end
    return findall(support)
end

"""
    refinecone(rel::EnvironmentRelations, cert; depth=rel.depth+1, optimizer=nothing,
               method=:auto, kwargs...)

The outer cone at `depth`, computed from a [`certifyrays`](@ref) certificate by wedge-restricted
discovery instead of the full deeper relations.

Every new extreme ray of the deeper cone lies outside the survivors' hull (an extreme ray that
decomposes over other cone members must be one of them), so only environments able to carry
weight beyond the survivor cone's violated facets can matter — together with the unwitnessed
survivors' fiber supports, which must re-emerge from the restricted system since their survival
is float-certified only. The witnessed rays seed the hull. With no eliminations the
certificate's cone is returned unchanged (`stable`). Remaining keyword arguments are passed to
`refinementrelations`.
"""
function refinecone(rel::EnvironmentRelations, cert; depth=rel.depth + 1, optimizer=nothing,
                    verbose::Bool=false, method::Symbol=:auto, kwargs...)
    _requireuniform(rel, "refinecone")
    eliminated = [Rational{BigInt}.(v.ray) for v in cert.verdicts if v.status === :eliminated]
    isempty(eliminated) && return cert.cone

    witnessed = [Rational{BigInt}.(v.ray) for v in cert.verdicts
                 if v.status === :finite || v.status === :periodic]
    undetermined = [Rational{BigInt}.(v.ray) for v in cert.verdicts
                    if v.status === :undetermined]
    d = size(rel.projection, 1)

    S = _tomatrix(vcat(witnessed, undetermined), d)
    length(_pivotcols!(copy(S))) == d ||
        throw(ArgumentError("the surviving rays do not span the composition space; " *
                            "compute the deeper cone directly"))

    t = @elapsed ((Ah, _, _) = facetsof(S; rays=true, incidence=false))
    verbose && (println("refinecone: ", length(eliminated), " eliminated rays, ",
                        size(Ah, 1), " survivor-hull facets (", round(t; digits=1), "s)");
                flush(stdout))
    keep = Set{Int}()
    ncut = 0
    for h in eachrow(Ah)
        any(r -> sum(h .* r) > 0, eliminated) || continue   # only facets the damage violates
        ncut += 1
        # survivor-hull facet normals carry huge integer entries; the exact slab LPs grind on
        # them, so with an optimizer the sweep runs in floats — over-collecting is sound
        t = @elapsed members = optimizer === nothing ? _wedgesupport(rel, collect(h)) :
                    _wedgemembers(optimizer, rel.relations, rel.projection, collect(h))
        union!(keep, members)
        verbose && (println("  cut facet ", ncut, ": ", length(members), " wedge members (",
                            round(t; digits=1), "s, keep now ", length(keep), ")");
                    flush(stdout))
    end
    for r in undetermined
        union!(keep, fibersupport(rel, r))
    end
    isempty(keep) &&
        throw(ArgumentError("empty wedge support; the certificate is inconsistent with `rel`"))
    verbose && (println("  keep: ", length(keep), " of ", length(rel.envs), " depth-",
                        rel.depth, " environments"); flush(stdout))

    t = @elapsed refined = refinementrelations(rel, sort!(collect(keep)); depth, kwargs...)
    verbose && (println("  refinement build: ", length(refined.envs), " environments, ",
                        size(refined.relations, 1), " relations (", round(t; digits=1), "s)");
                flush(stdout))
    t = @elapsed O = outercone(refined; method, optimizer, seeds=_tomatrix(witnessed, d))
    verbose && (println("  cone: ", size(O.rays, 1), " rays, ", size(O.facets, 1),
                        " facets (", round(t; digits=1), "s)"); flush(stdout))
    return O
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
    _requireuniform(rel, "refinementrelations")
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

_proportional(a, b) = !iszero(a) &&
                      all(a[i] * b[j] == a[j] * b[i]
                          for i in eachindex(a) for j in (i + 1):lastindex(a))

"""
    findraywitness(rel::EnvironmentRelations, m; maxscale=1, maxsize=Inf, periodic=false, nreps=2)

Search for a structure whose symmetrized composition is proportional to `m`, by reverse search
pruned to compositions elementwise at most `maxscale` times the primitive representative of `m`.
Returns the witness `Polyform`, or `nothing` if the pruned search is exhausted.

With `periodic=true`, structures are also accepted as unit cells: if some tiling of the
structure (see `Roly.tilings`; partial closures included, `nreps` passed through) brings the
*per-cell* composition — bulk bonds plus the bonds the tiling closes — onto `m`, the structure
is returned. Such a witness realizes `m` as an infinite periodic assembly.

A witness proves `m` realizable. Absence proves nothing — the cutoffs `maxsize`/`maxscale` are
inherent meta-parameters: a larger witness may always exist beyond them.
"""
function findraywitness(rel::EnvironmentRelations, m::AbstractVector; maxscale::Integer=1,
                        maxsize=Inf, periodic::Bool=false, nreps::Integer=2)
    d = size(rel.projection, 1)
    length(m) == d ||
        throw(DimensionMismatch("`m` has $(length(m)) entries but the composition space has $d"))
    mi = _primitive(m)
    bound = maxscale .* mi

    ns = nspecies(rel.rules)
    classmap = _bondclassmap(rel)
    hit = Ref{Any}(nothing)
    function f(s, _)
        c = _symcomposition(rel, classmap, s)
        all(c .<= bound) || return REJECT
        if _proportional(c, mi)
            hit[] = copy(s)
            return BREAK
        end
        if periodic
            for t in Roly.tilings(s; nreps)
                cc = copy(c)
                for β in t.bondtypes
                    cc[ns + classmap[β]] += 1
                end
                if _proportional(cc, mi)
                    hit[] = copy(s)
                    return BREAK
                end
            end
        end
        return ACCEPT
    end
    polyenum(f, rel.rules; maxsize)
    return hit[]
end

# One shared reverse search witnessing many rays at once: the pruning envelope is the union of
# the still-active rays' composition boxes and tightens as rays get witnessed, and each visited
# structure runs the tiling check once for all rays. A structure proportional to a ray witnesses
# it regardless of whose box admitted the visit — the boxes only steer pruning.
function _findraywitnesses(rel::EnvironmentRelations, rays; maxscale, maxsize, periodic, nreps)
    ns = nspecies(rel.rules)
    classmap = _bondclassmap(rel)
    mis = [_primitive(m) for m in rays]
    bounds = [maxscale .* mi for mi in mis]
    active = trues(length(rays))
    witnesses = Vector{Any}(nothing, length(rays))

    function f(s, _)
        c = _symcomposition(rel, classmap, s)
        any(i -> active[i] && all(c .<= bounds[i]), eachindex(bounds)) || return REJECT
        for i in eachindex(mis)
            active[i] || continue
            if _proportional(c, mis[i])
                witnesses[i] = copy(s)
                active[i] = false
            end
        end
        if periodic && any(active)
            for t in Roly.tilings(s; nreps)
                cc = copy(c)
                for β in t.bondtypes
                    cc[ns + classmap[β]] += 1
                end
                for i in eachindex(mis)
                    active[i] || continue
                    if _proportional(cc, mis[i])
                        witnesses[i] = copy(s)
                        active[i] = false
                    end
                end
            end
        end
        return any(active) ? ACCEPT : BREAK
    end
    polyenum(f, rel.rules; maxsize)
    return witnesses
end

_primitive(m::AbstractVector) = begin
    mq = Rational{BigInt}.(m)
    mi = numerator.(mq .* lcm(denominator.(mq)))
    mi .÷ reduce(gcd, mi)
end

"""
    certifyrays(rel::EnvironmentRelations; optimizer=nothing, maxscale=1, maxsize=6,
                periodic=true, refine=true, nreps=2, refinedepth=rel.depth+1, method=:auto)

Walk every extreme ray of the outer cone through the escalation ladder: finite witness,
periodic witness, fiber-restricted elimination one depth deeper.

  - returns `(; cone, verdicts, certified, stable)`: the outer cone, one
    `(ray, status, witness)` per extreme ray, whether `C = O_k` is certified, and whether the
    hierarchy has stabilized
  - `status` is `:finite` or `:periodic` (witnessed, with the witness structure or unit cell),
    `:eliminated` (provably not in the depth-`refinedepth` cone, so `O_k` is strictly loose
    there), or `:undetermined`
  - `certified == true` — every ray witnessed — proves the outer cone equals the true
    composition cone
  - `stable == true` — no ray eliminated — proves `O_refinedepth = O_k` without computing the
    deeper cone: witnessed rays lie in `C` and every other ray survived the refinement test, so
    the deeper cone contains all of this one's rays. Escalating depth is then pointless;
    only witness budget can settle the undetermined rays (`nothing` when `refine=false`)

This is a single pass at fixed parameters — nothing escalates automatically. The cone is
computed at `rel`'s depth and never changes; unwitnessed rays get exactly one elimination test
at `refinedepth`. `:undetermined` may flip to witnessed on a rerun with larger
`maxsize`/`maxscale`, or to `:eliminated` with a deeper `refinedepth` — the caller chooses
where to spend, per ray if need be, via the primitives. The elimination pass restricts one
shared refinement build to the union of the unwitnessed rays' fiber supports, which covers each
individual fiber and so stays lossless per ray.
"""
function certifyrays(rel::EnvironmentRelations; optimizer=nothing, maxscale::Integer=1,
                     maxsize=6, periodic::Bool=true, refine::Bool=true, nreps::Integer=2,
                     refinedepth::Integer=rel.depth + 1, seeds=nothing, method::Symbol=:auto)
    O = outercone(rel; method, optimizer, seeds)
    classmap = _bondclassmap(rel)

    rays = [collect(r) for r in eachrow(O.rays)]
    statuses = fill(:undetermined, length(rays))
    witnesses = _findraywitnesses(rel, rays; maxscale, maxsize, periodic, nreps)
    supports = Dict{Int,Vector{Int}}()
    for (i, m) in enumerate(rays)
        if witnesses[i] !== nothing
            finite = _proportional(_symcomposition(rel, classmap, witnesses[i]), _primitive(m))
            statuses[i] = finite ? :finite : :periodic
            continue
        end
        refine || continue
        U = fibersupport(rel, m; optimizer)
        # an empty fiber means the ray is not even in the exact depth-k cone (a float-mode
        # rationalization artifact); it is eliminated outright
        isempty(U) ? (statuses[i] = :eliminated) : (supports[i] = U)
    end

    if !isempty(supports)
        keep = sort!(union(values(supports)...))
        refined = refinementrelations(rel, keep; depth=refinedepth)
        for i in keys(supports)
            realizablecounts(refined, rays[i]; optimizer) === nothing &&
                (statuses[i] = :eliminated)
        end
    end

    verdicts = [(ray=rays[i], status=statuses[i], witness=witnesses[i],
                 support=get(supports, i, Int[])) for i in eachindex(rays)]
    certified = all(s -> s === :finite || s === :periodic, statuses)
    # every ray of O_k surviving at `refinedepth` (witnessed rays survive trivially: they are in
    # C) proves O_refinedepth = O_k — computing the deeper cone would return this same cone, so
    # undetermined rays are a witness-budget problem, not a depth problem
    stable = refine ? !any(s -> s === :eliminated, statuses) : nothing
    return (; cone=O, verdicts, certified, stable)
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
    intersect(a::OuterCone, b::OuterCone)

The intersection of two cones: an outer bound whenever both are, so cones from independent
(e.g. sampled) runs can be combined. Merging the underlying environment listings and rebuilding
is always at least as tight; intersection is the cheap route when only the cones are at hand.
"""
function Base.intersect(a::OuterCone, b::OuterCone)
    d = size(a.facets, 2)
    size(b.facets, 2) == d ||
        throw(DimensionMismatch("the cones live in different composition spaces"))
    A = [Rational{BigInt}.(a.facets); Rational{BigInt}.(b.facets)]
    A2, b2, _ = removeredundancy(A, zeros(Rational{BigInt}, size(A, 1)))
    rays, _ = extremerays(A2, b2)
    return _narrowcone(A2, rays)
end

"""
    outercone(rel::EnvironmentRelations; method=:auto, optimizer=nothing, seeds=nothing)

Compute the outer composition cone `O = Π(𝒦)` exactly, where `𝒦 = {μ ≥ 0 : relations ⋅ μ = 0}`
is the feasible cone of the relations.

  - `method`: `:rays` enumerates the extreme rays of `𝒦` and projects them — exact but hopeless
    once `𝒦` has many extreme rays; `:project` gift-wraps `O` directly with one LP per face
    certificate, so its cost scales with the size of `O`, not of `𝒦`; `:auto` picks by size
  - `optimizer`: a JuMP-compatible optimizer (e.g. `HiGHS.Optimizer`; requires JuMP loaded) to
    solve the `:project` LPs in floating point with rationalized results — much faster, but the
    cone is only as exact as the rationalization
  - `seeds`: known compositions in the cone, one per row, to start the `:project` hull from.
    The gift-wrap is output-sensitive, so seeding with the certified rays of a shallower cone
    concentrates the probes on the regions the deeper relations actually change

Soundness: the counts of every finite assembly lie in `𝒦`, so every composition lies in `O`.
Feasible points of `𝒦` are not claimed realizable; only the outer bound is exact.
"""
function outercone(rel::EnvironmentRelations; method::Symbol=:auto, optimizer=nothing,
                   seeds=nothing)
    ne = length(rel.envs)
    ns = size(rel.slacks, 1)
    nr = size(rel.relations, 1)
    if method === :auto
        method = optimizer === nothing && ne <= 64 ? :rays : :project
    end
    linearity = collect((ne + ns + 1):(ne + ns + nr))

    if method === :rays
        A = vcat(-Matrix{Rational{Int}}(I, ne, ne), Rational{Int}.(rel.slacks),
                 Rational{Int}.(rel.relations))
        rays, _ = extremerays(A, zeros(Rational{Int}, ne + ns + nr); linearity)
        projected = Rational{Int}.(rays) * transpose(rel.projection)
        facets, bf, _ = facetsof(projected; rays=true, incidence=false)
        minimalrays, _ = extremerays(Rational{Int}.(facets), Rational{Int}.(bf))
        return _narrowcone(facets, minimalrays)
    elseif method === :project
        # environments whose relation column vanishes are feasible on their own (the monomer
        # environments always are); their projections seed the hull for free
        # sparse: at large nₑ the dense identity block alone would be unallocatable
        A = [-sparse(one(Rational{Int}) * I, ne, ne); sparse(Rational{Int}.(rel.slacks));
             sparse(Rational{Int}.(rel.relations))]
        # environments with a vanishing relation column are feasible alone (monomers always
        # are); their projections seed the hull, joined by any caller-supplied compositions
        seedcols = findall(j -> nnz(rel.relations[:, j]) == 0, 1:ne)
        allseeds = permutedims(rel.projection[:, seedcols])
        if seeds !== nothing
            extra = Matrix{Rational{BigInt}}(undef, size(seeds, 1), size(allseeds, 2))
            for (i, r) in enumerate(eachrow(seeds))
                extra[i, :] .= Rational{BigInt}.(collect(r))
            end
            allseeds = vcat(Rational{BigInt}.(allseeds), extra)
        end
        facets, rays = projectcone(rel.projection, A; linearity, seeds=allseeds, optimizer)
        return _narrowcone(facets, rays)
    end
    throw(ArgumentError("unknown method `$(repr(method))`; expected `:auto`, `:rays` or `:project`"))
end
