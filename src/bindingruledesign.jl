### Design of the binding rules themselves, as opposed to `design.jl`, which takes the rules as given
### and designs the thermodynamic parameters. A design here is a `BindingRules`. The search starts
### from the fully addressable design, in which every site of the target carries its own color, and
### repeatedly merges pairs of species, which only ever adds structures.

"""
    fully_addressable_design(target::Polyform, particlespecies::AbstractVector{<:ParticleSpecies})
    fully_addressable_design(target::Polyform, particlespecies::ParticleSpecies)

Return the `BindingRules` in which every binding site of `target` is its own color, so the only
allowed bonds are the ones `target` actually makes.

This is the seed of the design search: it assembles `target` and nothing else, at the cost of using
as many colors as there are sites. Every particle of `target` becomes its own species, so
`particlespecies` is one species per particle, and passing a single `ParticleSpecies` uses it for all
of them.

Built from `bonds(target)` rather than from the adjacency matrix of `graphrep(target)`, because the
graph is stored in canonical vertex order, which need not follow particle order.
"""
function fully_addressable_design(target::Polyform, particlespecies::AbstractVector{<:ParticleSpecies})
    length(particlespecies) == nparticles(target) ||
        throw(ArgumentError("expected one particle species per particle of the target"))

    bondlist = Matrix{Int}(undef, nbonds(target), 4)
    for (k, bnd) in enumerate(bonds(target))
        lhs, rhs = first(bnd), last(bnd)
        bondlist[k, :] .= (lhs.particle, lhs.site, rhs.particle, rhs.site)
    end
    return BindingRules(bondlist, particlespecies)
end
fully_addressable_design(target::Polyform, particlespecies::ParticleSpecies) =
    fully_addressable_design(target, fill(particlespecies, nparticles(target)))

"""
    merge_species(rules::BindingRules, spcs_i, spcs_j, Δψ)

Return the rules obtained by identifying species `spcs_j` with species `spcs_i`, rotated by `Δψ`
sites.

Site `k` of `spcs_i` is merged with site `mod1(Δψ + k, n)` of `spcs_j`, where `n` is the number of
sites of `spcs_j`, and `spcs_j` is dropped. The merged site inherits every bond either site had, so
the merged rules allow at least every structure the original rules allowed.
"""
function merge_species(rules::BindingRules, spcs_i::Integer, spcs_j::Integer, Δψ::Integer)
    nsites_i = nsites(species(rules, spcs_i))
    nsites_j = nsites(species(rules, spcs_j))

    sites_i = 1:nsites_i
    sites_j = mod1.(Δψ .+ sites_i, nsites_j)

    colors_i = [color(rules, SpeciesSiteLoc(spcs_i, k)) for k in sites_i]
    colors_j = [color(rules, SpeciesSiteLoc(spcs_j, k)) for k in sites_j]

    imat = BitMatrix(interactionmatrix(rules))
    for (ci, cj) in zip(colors_i, colors_j)
        imat[ci, :] .|= imat[cj, :]
        imat[:, ci] .|= imat[:, cj]
    end

    keep = setdiff(axes(imat, 1), colors_j)
    imat = imat[keep, keep]

    geoms = collect(species(rules))
    deleteat!(geoms, spcs_j)
    return BindingRules(imat, geoms)
end

"""
    _canonicalforms(rules::BindingRules)

Return the canonical graph of `rules` together with that of its mirror image.

Two rule sets allow the same structures, that is they agree up to relabeling and rotating the
species, exactly when their canonical graphs are equal, so comparing these forms decides equivalence
outright rather than through an id. The mirror is the transposed graph, since mirroring reverses the
site cycle of every species. Storing both forms but querying only the first lets a single lookup
recognize a rule set already seen in either handedness.
"""
function _canonicalforms(rules::BindingRules)
    g = graphrep(rules)
    mirror = NautyDiGraph(permutedims(adjacency_matrix(g)); vertex_labels=labels(g))
    canonize!(g)
    canonize!(mirror)
    return g, mirror
end

"""
    generate_refinements(rules::BindingRules; enable_rotations=true)

Return every rule set obtained by merging one pair of species of `rules` that carry the same number
of sites, over all relative rotations unless `enable_rotations` is false.

Rotations may only be disabled when every particle carries the same orientation, since otherwise the
zero-rotation merge is not the canonical one.
"""
function generate_refinements(rules::BindingRules; enable_rotations=true)
    m = nspecies(rules)
    offspring = typeof(rules)[]
    m == 1 && return offspring

    for i in 1:m, j in i+1:m
        nsites_i = nsites(species(rules, i))
        nsites_j = nsites(species(rules, j))
        nsites_i == nsites_j || continue

        if enable_rotations
            for Δψ in 0:nsites_i-1
                push!(offspring, merge_species(rules, i, j, Δψ))
            end
        else
            push!(offspring, merge_species(rules, i, j, 0))
        end
    end
    return offspring
end

"""
    shapeid(p::Polyform)

Return the canonical id of `p` with its species labels erased, so that two polyforms share an id
exactly when they have the same shape regardless of which species sit where.
"""
function shapeid(p::Polyform)
    g = copy(graphrep(p))
    setlabels!(g, zeros(Cint, nv(g)))
    return canonical_id(g)
end

"""
    find_shape(p::Polyform, strs::AbstractVector{<:Polyform})

Return the positions in `strs` of the structures that have the same shape as `p`, ignoring species
assignment.

`strs` is assumed sorted by number of particles, as `polygen` returns it, so the scan stops once the
structures grow past `p`.
"""
function find_shape(p::Polyform, strs::AbstractVector{<:Polyform})
    needle = shapeid(p)
    n = nparticles(p)

    matches = Int[]
    for (i, str) in enumerate(strs)
        m = nparticles(str)
        m > n && break
        m < n && continue
        shapeid(str) == needle && push!(matches, i)
    end
    return matches
end

"""
    enumerate_structures(rules::BindingRules; maxsize, maxpolys)

Return `(strs, valid)`, the structures `rules` allows and whether the enumeration closed strictly
inside both bounds.

Hitting either bound means the rule set is unbounded as far as the enumeration can tell, and such a
rule set is rejected rather than scored.
"""
function enumerate_structures(rules::BindingRules; maxsize, maxpolys)
    strs = polygen(rules; maxsize, maxpolys)
    valid = length(strs) < maxpolys && maximum(nparticles, strs) < maxsize
    return strs, valid
end

"""
    evaluate_system(rules::BindingRules, target::Polyform, ε::Real, ϕ::Real; maxsize, maxpolys)

Return `(valid, score)`, scoring `rules` by the equilibrium yield of `target` at bond energy `ε` and
total particle concentration `ϕ`, with every bond type held at the same energy.

`valid` is false, and the score zero, if the rule set is unbounded or cannot build `target`. All
structures are given unit partition function, so the score depends on the rules alone and no entropy
model is needed.
"""
function evaluate_system(rules::BindingRules, target::Polyform, ε::Real, ϕ::Real; maxsize, maxpolys)
    strs, valid = enumerate_structures(rules; maxsize, maxpolys)
    valid || return false, 0.0

    target_idxs = find_shape(target, strs)
    isempty(target_idxs) && return false, 0.0

    ns = nspecies(rules)
    M = permutedims(reduce(hcat, composition.(strs)))
    # One uniform bond energy, so the bond columns collapse into a single total-bond column.
    M = hcat(M[:, 1:ns], sum(M[:, ns+1:end]; dims=2))

    Ωs = ones(size(M, 1))
    target_yields = map(target_idxs) do i
        ϕs = ϕ / nparticles(target) * M[i, 1:end-1]
        yields(ϕs, [ε], M, Ωs)[i]
    end
    return true, maximum(target_yields)
end

"""
    evaluate_system(rules::BindingRules; maxsize, maxpolys, use_enum=false, no_opposing_particles=false, target=nothing)

Return `(valid, score)`, scoring `rules` by how few structures it allows, so that the search prefers
the most selective rule set.

Any rule set whose enumeration runs past `maxsize` or `maxpolys` is rejected, as is one that cannot
build `target`, which may be a single `Polyform` or a collection of them. Set `use_enum` to count
structures without materializing them, which is faster but then `target` cannot be checked.
"""
function evaluate_system(rules::BindingRules; maxsize, maxpolys, use_enum=false,
                         no_opposing_particles=false, target=nothing)
    if use_enum
        result = polyenum(rules; maxsize, maxpolys)
        nstrs = result.n
        valid = nstrs < maxpolys && result.maxsize < maxsize
        strs = nothing
    else
        strs, valid = enumerate_structures(rules; maxsize, maxpolys)
        nstrs = length(strs)
    end

    if valid && !use_enum && !isnothing(target)
        targets = target isa AbstractArray ? target : (target,)
        target_idxs = Int[]
        for t in targets
            idxs = find_shape(t, strs)
            isempty(idxs) && return false, 0.0
            append!(target_idxs, idxs)
        end

        if no_opposing_particles
            valid = any(!has_opposing_particles(strs[t]) for t in target_idxs)
        end
    end

    return valid ? (true, convert(Float64, -nstrs)) : (false, 0.0)
end

"""
    has_opposing_particles(str::Polyform)

Return whether `str` contains particles at more than two distinct orientations, or at two that
differ by half a turn.

Only meaningful for squares, and reaches into the particle poses directly.
"""
function has_opposing_particles(str::Polyform)
    # TODO: generalize beyond squares and off the internal pose representation.
    angles = (rotation_angle(part.pose.psi) for part in str.particles)
    orientations = sort!(unique(x -> round(x; sigdigits=3), collect(angles)))
    nos = length(orientations)
    return nos > 2 || (nos == 2 && orientations[2] == orientations[1] + 1)
end

"""
    _search_designs(seed::BindingRules, evaluator; greed, enable_rotations, verbose)

Return `(systems, network)` from a breadth-first search over the refinements of `seed`, keeping those
`evaluator` accepts and, when `greed` is given, only its `greed` highest scoring children at each
step.

`evaluator` maps a `BindingRules` to `(valid, score)`. `network` records which rule set was refined
from which, as a graph on `1:length(systems)`.
"""
function _search_designs(seed::BindingRules, evaluator; greed, enable_rotations, verbose)
    systems = typeof(seed)[seed]
    seen = Set(_canonicalforms(seed))
    system_network = Graph(1)

    # Plain vector as a FIFO: entries are only appended and consumed in order.
    queue = Tuple{Int,typeof(seed)}[(1, seed)]
    head = 1

    while head <= length(queue)
        nsys = length(systems)
        n, sys = queue[head]
        head += 1
        if verbose
            print("Evaluating system[$(nspecies(sys)), $(nbonds(sys))], nsys = $nsys, ",
                  "$(length(queue) - head + 1) in queue            \r")
        end

        new_systems = generate_refinements(sys; enable_rotations)
        new_forms = _canonicalforms.(new_systems)
        unique_idxs = unique(i -> first(new_forms[i]), eachindex(new_forms))
        new_systems, new_forms = new_systems[unique_idxs], new_forms[unique_idxs]

        fresh = [first(f) ∉ seen for f in new_forms]
        new_systems, new_forms = new_systems[fresh], new_forms[fresh]

        scores = zeros(length(new_systems))
        valids = zeros(Bool, length(new_systems))
        @threads for i in eachindex(new_systems)
            valids[i], scores[i] = evaluator(new_systems[i])
        end
        new_systems, new_forms, scores = new_systems[valids], new_forms[valids], scores[valids]

        if !isnothing(greed)
            ranking = sortperm(scores; rev=true)
            length(ranking) > greed && (ranking = ranking[1:greed])
            new_systems, new_forms = new_systems[ranking], new_forms[ranking]
        end

        for (j, (newsys, newform)) in enumerate(zip(new_systems, new_forms))
            push!(systems, newsys)
            push!(seen, newform...)
            push!(queue, (nsys + j, newsys))

            add_vertex!(system_network)
            add_edge!(system_network, n, nsys + j)
        end
    end

    verbose && println("\n\nGeneration complete: generated $(length(systems)) systems.")
    return systems, system_network
end

"""
    create_designs(target::Polyform, particlespecies; maxsize, maxpolys, greed=nothing, ε=nothing, ϕ=nothing, enable_rotations=true, verbose=false, kwargs...)

Return `(systems, network)`, the rule sets that assemble `target`, found by starting from the fully
addressable design and merging species as far as the allowed structures stay bounded.

Supplying both `ε`, a uniform bond energy, and `ϕ`, a total particle concentration, scores designs by
the equilibrium yield of the target rather than by the number of structures they allow. `greed` keeps
only that many highest scoring children per step, and `maxsize` and `maxpolys` bound the enumeration
used to decide whether a rule set is still bounded.
"""
function create_designs(target::Polyform, particlespecies;
                        maxsize::Integer,
                        maxpolys::Integer,
                        greed::Union{Integer,Nothing}=nothing,
                        ε::Union{Real,Nothing}=nothing,
                        ϕ::Union{Real,Nothing}=nothing,
                        enable_rotations=true,
                        verbose=false,
                        kwargs...)
    # TODO: all particles need the same orientation for rotations to be disabled safely.
    seed = fully_addressable_design(target, particlespecies)

    evaluator = if !isnothing(ε) && !isnothing(ϕ)
        sys -> evaluate_system(sys, target, ε, ϕ; maxsize, maxpolys, kwargs...)
    else
        sys -> evaluate_system(sys; maxsize, maxpolys, target, kwargs...)
    end

    return _search_designs(seed, evaluator; greed, enable_rotations, verbose)
end

"""
    create_multifarious_designs(rules::BindingRules, targets; maxsize, maxpolys, greed=nothing, enable_rotations=true, verbose=false, kwargs...)

Return `(systems, network)` as [`create_designs`](@ref), but starting from an existing rule set and
requiring every structure in `targets` to remain assemblable.
"""
function create_multifarious_designs(rules::BindingRules, targets;
                                     maxsize::Integer,
                                     maxpolys::Integer,
                                     greed::Union{Integer,Nothing}=nothing,
                                     enable_rotations=true,
                                     verbose=false,
                                     kwargs...)
    evaluator = sys -> evaluate_system(sys; maxsize, maxpolys, target=targets, kwargs...)
    return _search_designs(rules, evaluator; greed, enable_rotations, verbose)
end
