"""
    EnvironmentRelations

The exact linear relations satisfied by the particle-environment counts of every finite assembly,
at a fixed environment radius.

  - `envs`: every particle environment the rules allow; the count vector `μ` is indexed by it
  - `relations`: `relations * μ == 0` exactly for the counts of any finite assembly or soup
  - `projection`: `projection * μ == [n_α; b̃]`, the composition with bond counts grouped by
    symmetry-equivalent bond type
  - `bondclasses`: the `((species, label), (species, label))` pair labelling each `b̃` axis
  - `depth`: environment radius

Each undirected bond is seen twice, once from each endpoint, and both sightings crop to the same
bond neighborhood with the endpoint order swapped (`Roly.rootbonds`), so the sightings of a crop
and of its reversal count the same bonds. This gives one exact integer relation per asymmetric
crop pair, with the bond count itself eliminated: the slot system `Σ_w N_e(w) μ_w = z_c =
Σ_w N_ē(w) μ_w` becomes `Σ_w (N_e(w) - N_ē(w)) μ_w = 0`, and `z_c ≥ 0` is implied. A crop equal
to its own reversal balances itself and contributes no row. All inequality content of the cone
`{μ ≥ 0 : relations * μ = 0}` lives in `μ ≥ 0`.

The equalities are exact for the unsaturated ensemble, where unbound sites are allowed
everywhere. A saturated ensemble would need the same equalities and a different admissibility
convention for the environments; it is not implemented.
"""
struct EnvironmentRelations{E,S}
    envs::Vector{E}
    relations::SparseMatrixCSC{Int,Int}
    projection::Matrix{Rational{Int}}
    bondclasses::Vector{NTuple{2,Tuple{Int,Int}}}
    depth::Int
    rules::S
end

Base.show(io::IO, rel::EnvironmentRelations) =
    print(io, "EnvironmentRelations[k=", rel.depth, ", ", length(rel.envs), " environments, ",
          size(rel.relations, 1), " relations, d=", size(rel.projection, 1), "]")

_speciesranges(rules) = [extrema(labels(graphrep(Roly.species(rules, i)))) for i in 1:nspecies(rules)]
_speciesoflabel(ranges, l) = findfirst(r -> r[1] <= l <= r[2], ranges)
# Bond types up to particle symmetry: sites with equal labels are equivalent by Roly's label
# semantics, so the unordered pair of (species, label) is the finest bond resolution that is a
# class function of environments.
_siteclass(ranges, l) = (_speciesoflabel(ranges, l), l)

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
    isempty(envs) && throw(ArgumentError("no environments to build relations from"))
    offset = Roly._markoffset(rules)
    ranges = _speciesranges(rules)
    ns = nspecies(rules)
    ne = length(envs)

    bondclasses = prescribed === nothing ? NTuple{2,Tuple{Int,Int}}[] : collect(prescribed)
    clsidx = Dict(cls => i for (i, cls) in enumerate(bondclasses))
    rows = Dict{typeof(first(envs).graph),Int}()   # crop ⇒ row, keyed by the orientation seen first
    I, J, V = Int[], Int[], Int[]
    rootspecies = zeros(Int, ne)
    halfbonds = Tuple{Int,Int}[]                   # (bond class, environment) per root bond
    for (j, w) in enumerate(envs)
        rootspecies[j] = _speciesoflabel(ranges, mod1(labels(w.graph)[w.rootvertex], offset))
        for b in rootbonds(w, depth)
            cls = minmax(_siteclass(ranges, b.sites[1]), _siteclass(ranges, b.sites[2]))
            c = get!(clsidx, cls) do
                prescribed === nothing ||
                    throw(ArgumentError("the environments carry a bond class the prescribed axes lack"))
                push!(bondclasses, cls)
                length(bondclasses)
            end
            push!(halfbonds, (c, j))

            b.view == b.reversed && continue
            r, s = get(rows, b.view, 0), 1
            r == 0 && (r = get(rows, b.reversed, 0); s = -1)
            r == 0 && (r = rows[b.view] = length(rows) + 1; s = 1)
            push!(I, r)
            push!(J, j)
            push!(V, s)
        end
    end
    # `sparse` sums repeated entries, so this is N_e(w) - N_ē(w) per row and column
    relations = dropzeros!(sparse(I, J, V, length(rows), ne))

    # order the bond axes deterministically; prescribed axes keep their given order
    order = prescribed === nothing ? sortperm(bondclasses) : collect(eachindex(bondclasses))
    axis = invperm(order)
    projection = zeros(Rational{Int}, ns + length(bondclasses), ne)
    for j in 1:ne
        projection[rootspecies[j], j] = 1
    end
    for (c, j) in halfbonds
        projection[ns + axis[c], j] += 1//2
    end
    return EnvironmentRelations(envs, relations, projection, bondclasses[order], Int(depth), rules)
end

"""
    environmentcounts(rel::EnvironmentRelations, poly::Polyform)

Count the particles of `poly` into `rel`'s environment bins, as a vector `μ` indexed like
`rel.envs`. Every particle must match some bin, which holds whenever the relations were built
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

# The fiber `{μ ≥ 0 : relations⋅μ = 0, projection⋅μ = m}` as an H-representation `A*x <= b` with
# linearity rows. It is a polytope: every environment contributes +1 to its root-species row, so
# `Πμ = m` bounds `Σμ`, hence all LPs over it are bounded.
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
    realizablecounts(rel::EnvironmentRelations, m; optimizer=nothing)

A count vector `μ ≥ 0` with `relations⋅μ = 0` and `projection⋅μ = m`, or `nothing` when none
exists. Feasibility certifies nothing about assembly; `nothing` rules `m` out exactly.

Without an `optimizer` the fiber LP is solved exactly (requires Polymake loaded), which is viable
only for small systems. With one (requires JuMP loaded), the LP is solved in floating point and an
infeasibility verdict is only reported after its Farkas certificate verifies in exact
arithmetic; an unverifiable certificate returns `missing` (unknown) instead, so `nothing`
always remains an exact elimination.
"""
function realizablecounts(rel::EnvironmentRelations, m::AbstractVector; optimizer=nothing)
    if optimizer !== nothing
        d = size(rel.projection, 1)
        length(m) == d ||
            throw(DimensionMismatch("`m` has $(length(m)) entries but the composition space has $d"))
        return _fiberfeasible(optimizer, rel.relations, rel.projection, m)
    end
    length(rel.envs) <= 4096 ||
        throw(ArgumentError("the exact fiber LP is infeasible at $(length(rel.envs)) " *
                            "environments; pass a JuMP `optimizer`"))
    A, b, linearity = _fibersystem(rel, m)
    status, x, _ = solvelp(A, b, zeros(Rational{Int}, length(rel.envs)); linearity)
    return status === :optimal ? x : nothing
end

_fiberfeasible(optimizer, L, P, m) =
    throw(ArgumentError("passing an `optimizer` requires JuMP to be loaded, e.g. `using JuMP, HiGHS`"))

_fibermembers(optimizer, L, P, m) =
    throw(ArgumentError("passing an `optimizer` requires JuMP to be loaded, e.g. `using JuMP, HiGHS`"))

"""
    fibersupport(rel::EnvironmentRelations, m; optimizer=nothing)

Indices of every environment that can carry weight in a count vector realizing the composition
`m`: the support of the fiber `{μ ≥ 0 : relations⋅μ = 0, projection⋅μ = m}`. Empty when `m` is
infeasible.

Cropping maps any deeper fiber over `m` into this one, so only refinements of these environments
can matter when testing whether `m` survives at a larger depth; see [`refinementrelations`](@ref).

Without an `optimizer` the scan runs one exact LP per environment (requires Polymake loaded),
which dominates certification on any system of size. With one (requires JuMP loaded) the same
scan runs in floating point on a single model with the objective swapped per environment.
Over-collecting is sound here, since the support only decides which environments get refined, so
float membership needs no certificate.
"""
function fibersupport(rel::EnvironmentRelations, m::AbstractVector; optimizer=nothing)
    optimizer === nothing || return _fibermembers(optimizer, rel.relations, rel.projection, m)
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
        bindingsite(Roly.species(rules, loc[1]), loc[2]).vertices)]
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

# Whether a tiling of a structure with composition `c` can bring its per-cell composition onto the
# ray `mi`: a tiling adds bonds but no particles, so the species counts must already be
# proportional to the ray's, and every bond class must still be at or below its target count.
function _cellcandidate(c, mi, ns)
    p = findfirst(>(0), view(mi, 1:ns))
    p === nothing && return false
    c[p] > 0 || return false
    all(c[i] * mi[p] == c[p] * mi[i] for i in 1:ns) || return false
    return all(c[i] * mi[p] <= c[p] * mi[i] for i in (ns + 1):length(mi))
end

"""
    findraywitness(rel::EnvironmentRelations, m; maxscale=1, maxsize=Inf, periodic=false, nreps=2)

Search for a structure whose symmetrized composition is proportional to `m`, by reverse search
pruned to compositions elementwise at most `maxscale` times the primitive representative of `m`.
Returns the witness `Polyform`, or `nothing` if the pruned search is exhausted.

With `periodic=true`, structures are also accepted as unit cells: if some tiling of the
structure (see `Roly.tilings`; partial closures included, `nreps` passed through) brings the
per-cell composition, bulk bonds plus the bonds the tiling closes, onto `m`, the structure is
returned. Such a witness realizes `m` as an infinite periodic assembly. Tilings are only
enumerated for structures whose species counts are proportional to `m` and whose bond counts do
not exceed the matching multiple of `m`, since a tiling adds bonds but no particles.

A witness proves `m` realizable. Absence proves nothing: the cutoffs `maxsize` and `maxscale` are
inherent meta-parameters, and a larger witness may always exist beyond them.
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
        if periodic && _cellcandidate(c, mi, ns)
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
# it regardless of whose box admitted the visit; the boxes only steer pruning.
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
        if periodic && any(i -> active[i] && _cellcandidate(c, mis[i], ns), eachindex(mis))
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
    `(ray, status, witness, support)` per extreme ray, whether `C = O_k` is certified, and
    whether the hierarchy has stabilized
  - `status` is `:finite` or `:periodic` (witnessed, with the witness structure or unit cell),
    `:eliminated` (provably not in the depth-`refinedepth` cone, so `O_k` is strictly loose
    there), or `:undetermined`
  - `certified == true`, every ray witnessed, proves the outer cone equals the true composition
    cone
  - `stable == true`, no ray eliminated, proves `O_refinedepth = O_k` without computing the
    deeper cone: witnessed rays lie in `C` and every other ray survived the refinement test, so
    the deeper cone contains all of this one's rays. Escalating depth is then pointless;
    only witness budget can settle the undetermined rays (`nothing` when `refine=false`)

This is a single pass at fixed parameters; nothing escalates automatically. The cone is
computed at `rel`'s depth and never changes; unwitnessed rays get exactly one elimination test
at `refinedepth`. `:undetermined` may flip to witnessed on a rerun with larger
`maxsize`/`maxscale`, or to `:eliminated` with a deeper `refinedepth`. The elimination pass
restricts one shared refinement build to the union of the unwitnessed rays' fiber supports,
which covers each individual fiber and so stays lossless per ray.
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
    # C) proves O_refinedepth = O_k, so undetermined rays are a witness-budget problem, not a
    # depth problem
    stable = refine ? !any(s -> s === :eliminated, statuses) : nothing
    return (; cone=O, verdicts, certified, stable)
end

"""
    OuterCone

Outer bound on the composition cone: every composition of a finite assembly or soup satisfies
`facets * m ≤ 0`. Membership certifies nothing about realizability; the bound is one-sided.

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

# The rows of `seeds` that lie in the outer cone of `rel`, as a matrix. A seed proven outside throws, since
# the caller claimed it a member; a seed whose infeasibility certificate does not verify is dropped.
function _checkedseeds(rel::EnvironmentRelations, seeds, optimizer)
    d = size(rel.projection, 1)
    kept = Vector{Rational{BigInt}}[]
    for (i, r) in enumerate(eachrow(seeds))
        s = Rational{BigInt}.(collect(r))
        length(s) == d ||
            throw(DimensionMismatch("seed $i has $(length(s)) entries but the composition space has $d"))
        x = realizablecounts(rel, s; optimizer)
        if x === nothing
            throw(ArgumentError("seed $i ($(s)) lies outside the depth-$(rel.depth) outer cone, so it is " *
                                "not a realizable composition; check how it was certified"))
        elseif x === missing
            @warn "membership of seed $i in the outer cone could not be decided; dropping it"
        else
            push!(kept, s)
        end
    end
    return isempty(kept) ? zeros(Rational{BigInt}, 0, d) : permutedims(reduce(hcat, kept))
end

"""
    outercone(rel::EnvironmentRelations; method=:auto, optimizer=nothing, seeds=nothing,
              verbose=false)

Compute the outer composition cone `O = Π(𝒦)` exactly, where `𝒦 = {μ ≥ 0 : relations ⋅ μ = 0}`
is the feasible cone of the relations.

  - `method`: `:rays` enumerates the extreme rays of `𝒦` and projects them, exact and the faster
    route on small systems, but the ray count is exponential in the number of environments and
    the enumeration will consume all available memory before it gives up, so raise it past the
    `:auto` bound only on a system you have measured; `:project` gift-wraps `O` directly with one
    LP per face certificate, so its cost scales with the size of `O`, not of `𝒦`; `:auto` picks
    by size
  - `optimizer`: a JuMP-compatible optimizer (e.g. `HiGHS.Optimizer`; requires JuMP loaded) to
    solve the `:project` LPs in floating point with rationalized results; much faster, but the
    cone is only as exact as the rationalization
  - `seeds`: compositions known to lie in this cone, one per row, to start the `:project` hull
    from. The gift-wrap is output-sensitive, so seeding with the witnessed rays of a shallower
    cone concentrates the probes on the regions the deeper relations actually change. The hull
    only ever grows, so a seed outside the cone would enlarge the result; every seed is therefore
    checked with one fiber LP (see [`realizablecounts`](@ref)) before use. A seed proven outside
    the cone throws an `ArgumentError`, and a seed whose membership cannot be decided is dropped
    with a warning
  - `verbose`: report every support LP of the `:project` gift-wrap as it happens

Soundness: the counts of every finite assembly lie in `𝒦`, so every composition lies in `O`.
Feasible points of `𝒦` are not claimed realizable; only the outer bound is exact.
"""
function outercone(rel::EnvironmentRelations; method::Symbol=:auto, optimizer=nothing,
                   seeds=nothing, verbose::Bool=false)
    ne = length(rel.envs)
    nr = size(rel.relations, 1)
    if method === :auto
        # `:rays` enumerates the extreme rays of `𝒦`, which is exponential in nₑ and unbounded in
        # memory (raising the bound to 128 sent the enumeration to 100 GB on a mid-sized
        # fixture), and nₑ is a poor predictor of that ray count; callers who know their system
        # pass `method` explicitly
        method = ne <= 64 ? :rays : :project
    end
    linearity = collect((ne + 1):(ne + nr))

    if method === :rays
        A = vcat(-Matrix{Rational{Int}}(I, ne, ne), Rational{Int}.(rel.relations))
        rays, _ = extremerays(A, zeros(Rational{Int}, ne + nr); linearity)
        projected = Rational{Int}.(rays) * transpose(rel.projection)
        facets, bf, _ = facetsof(projected; rays=true, incidence=false)
        minimalrays, _ = extremerays(Rational{Int}.(facets), Rational{Int}.(bf))
        return _narrowcone(facets, minimalrays)
    elseif method === :project
        # sparse: at large nₑ the dense identity block alone would be unallocatable
        A = [-sparse(one(Rational{Int}) * I, ne, ne); sparse(Rational{Int}.(rel.relations))]
        # environments with a vanishing relation column are feasible alone (monomers always
        # are); their projections seed the hull, joined by any caller-supplied compositions
        seedcols = findall(j -> nnz(rel.relations[:, j]) == 0, 1:ne)
        allseeds = permutedims(rel.projection[:, seedcols])
        if seeds !== nothing
            extra = _checkedseeds(rel, seeds, optimizer)
            allseeds = vcat(Rational{BigInt}.(allseeds), extra)
        end
        facets, rays = projectcone(rel.projection, A; linearity, seeds=allseeds, optimizer, verbose)
        return _narrowcone(facets, rays)
    end
    throw(ArgumentError("unknown method `$(repr(method))`; expected `:auto`, `:rays` or `:project`"))
end
