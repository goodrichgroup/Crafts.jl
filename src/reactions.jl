"""
    Reaction

A single aggregation/fragmentation reaction `i + j ⟷ k`, as handed to a rate kernel. 

A lightweightview into a [`ReactionNetwork`](@ref); use [`indices`](@ref), [`reactants`](@ref), [`product`](@ref),
[`bondcounts`](@ref), [`cut`](@ref) and [`halves`](@ref) to read it.
"""
struct Reaction{N}
    network::N
    index::Int
end

"""
    indices(rxn::Reaction)

`(i, j, k)`: the positions of the two reactants and the product in the structure list.
"""
indices(rxn::Reaction) = rxn.network._reactions[rxn.index]

"""
    product(rxn::Reaction)

The structure that forms, carrying the geometry its fragments have once bound.
"""
product(rxn::Reaction) = structures(rxn.network.assemblysystem)[indices(rxn)[3]]

"""
    reactants(rxn::Reaction)

The two fragments, as the canonical structures of their isomorphism class. These describe the same
fragments as `product` restricted to [`halves`](@ref), but in their free rather than bound poses.
"""
function reactants(rxn::Reaction)
    strs = structures(rxn.network.assemblysystem)
    i, j, _ = indices(rxn)
    return (strs[i], strs[j])
end

"""
    bondcounts(rxn::Reaction)

How many bonds of each type are formed, equivalently broken, as a vector of length `nbonds`. Obtained
from the composition of the product less those of its fragments; the particle counts cancel, since
particles are conserved.
"""
function bondcounts(rxn::Reaction)
    asys = rxn.network.assemblysystem
    M = compositionmatrix(asys)
    i, j, k = indices(rxn)
    b = nspecies(asys)+1:size(M, 2)
    return M[k, b] - M[i, b] - M[j, b]
end

"""
    cut(rxn::Reaction)

The edges broken to separate the product into its two fragments.
"""
cut(rxn::Reaction) = rxn.network._cuts[rxn.index]

"""
    halves(rxn::Reaction)

The vertices of the product belonging to each fragment.
"""
halves(rxn::Reaction) = rxn.network._halves[rxn.index]

function Base.show(io::Core.IO, rxn::Reaction)
    i, j, k = indices(rxn)
    return print(io, "Reaction[$i + $j ⟷ $k, $(length(cut(rxn))) bond(s)]")
end

# An exterior edge is one whose reverse is also present; `src < dst` keeps one edge per bond.
_exterior_edges(g::AbstractGraph) = (e for e in edges(g) if has_edge(g, e.dst, e.src) && e.src < e.dst)

function cleave(g::AbstractGraph, es)
    g = copy(g)
    for edge in es
        revedge = reverse(edge)
        revedge ∈ edges(g) || throw(ArgumentError("Cannot cleave a nondirected (interior) edge."))
        rem_edge!(g, edge)
        rem_edge!(g, revedge)
    end
    return g, connected_components(g)
end

are_separated(vi, vj, vs1, vs2) = (vi ∈ vs1 && vj ∈ vs2) || (vi ∈ vs2 && vj ∈ vs1)

"""
    generate_cuts(g; maxbonds)

Every way of splitting `g` into exactly two connected pieces by breaking at most `maxbonds` bonds.
Returns the broken edges and the vertices of each piece.
"""
function generate_cuts(g::AbstractGraph; maxbonds)
    edgelist = collect(_exterior_edges(g))
    cuts = Vector{eltype(edgelist)}[]
    halves = Tuple{Vector{eltype(g)},Vector{eltype(g)}}[]
    isempty(edgelist) && return cuts, halves

    current_cut = [first(edgelist)]

    # Depth-first backtracking search over cuts of length <= maxbonds
    while true
        current_edge = current_cut[end]
        _, components = cleave(g, current_cut)

        if length(components) > 1
            # A reaction has exactly two reactants, so ignore cuts that shatter `g` further.
            if length(components) == 2 &&
               all(e -> are_separated(e.src, e.dst, components...), current_cut)
                push!(cuts, copy(current_cut))
                push!(halves, (components[1], components[2]))
            end
            nextedge_idx = nothing                                  # traverse back up
        elseif length(current_cut) == maxbonds
            nextedge_idx = nothing                                  # traverse back up
        else
            nextedge_idx = findfirst(>(current_edge), edgelist)     # keep going down
        end

        while isnothing(nextedge_idx)
            isempty(current_cut) && @goto finished
            current_edge = pop!(current_cut)
            nextedge_idx = findfirst(>(current_edge), edgelist)
        end
        push!(current_cut, edgelist[nextedge_idx])
    end
    @label finished

    return cuts, halves
end

"""
    ReactionNetwork(asys::AssemblySystem; maxbonds=Inf, fwdkernel=Returns(1.0), bwdkernel=fwdkernel)

The reactions `i + j ⟷ k` among the structures of `asys` that break at most `maxbonds` bonds at a
time, with a forward and a backward rate for each.

Each kernel takes a [`Reaction`](@ref) and returns a rate, read back with [`fwdrates`](@ref) and
[`bwdrates`](@ref). `active` marks the reactions that ended up with a nonzero rate; the others cannot
contribute and may be skipped.

Iterating a `ReactionNetwork` yields one `Reaction` per entry. Use [`rate!`](@ref) to apply different
kernels without repeating the cut enumeration, which is the expensive part.
"""
struct ReactionNetwork{A,E,V}
    assemblysystem::A
    _reactions::Vector{NTuple{3,Int}}
    _cuts::Vector{Vector{E}}
    _halves::Vector{Tuple{V,V}}

    _fwdrates::Vector{Float64}
    _bwdrates::Vector{Float64}
    active::Vector{Bool}
end

"""
    assemblysystem(net::ReactionNetwork)

The [`AssemblySystem`](@ref) whose structures `net` reacts.
"""
assemblysystem(net::ReactionNetwork) = net.assemblysystem

"""
    reactions(net::ReactionNetwork)
    fwdrates(net::ReactionNetwork)
    bwdrates(net::ReactionNetwork)

The `(i, j, k)` index triples of every reaction, and their forward and backward rates. Positions line
up with `net[r]`, with iterating `net`, and with each other.

These are the network's own storage rather than copies, so they allocate nothing but should be treated
as read-only; use [`rate!`](@ref) to change rates.
"""
reactions(net::ReactionNetwork) = net._reactions
fwdrates(net::ReactionNetwork) = net._fwdrates
bwdrates(net::ReactionNetwork) = net._bwdrates

"""
    isdetailedbalanced(net::ReactionNetwork)

Whether every active reaction has equal forward and backward rates.
"""
function isdetailedbalanced(net::ReactionNetwork)
    kfwd, kbwd = fwdrates(net), bwdrates(net)
    for r in eachindex(reactions(net))
        net.active[r] || continue
        kfwd[r] == kbwd[r] || return false
    end
    return true
end

# Fresh copies of the active subset, for consumers that rescale the rates in place.
_copy_active_reactions(net::ReactionNetwork) = net._reactions[net.active]
_copy_active_fwdrates(net::ReactionNetwork) = net._fwdrates[net.active]
_copy_active_bwdrates(net::ReactionNetwork) = net._bwdrates[net.active]

function ReactionNetwork(asys::AssemblySystem; maxbonds=Inf, kwargs...)
    gs = [graphrep(s) for s in structures(asys)] # assume all graphs are already canonized (which is true for Roly)
    ids = Dict(hash(g) => i for (i, g) in enumerate(gs))

    E = edgetype(first(gs))
    V = Vector{eltype(first(gs))}
    reactions = NTuple{3,Int}[]
    cuts = Vector{E}[]
    halves = Tuple{V,V}[]

    for (k, g) in enumerate(gs)
        _cuts, _halves = generate_cuts(g; maxbonds)
        for (c, half) in zip(_cuts, _halves)
            g1, g2 = g[half[1]], g[half[2]]
            canonize!(g1)
            canonize!(g2)
            push!(reactions, (ids[hash(g1)], ids[hash(g2)], k))
            push!(cuts, c)
            push!(halves, half)
        end
    end

    n = length(reactions)
    net = ReactionNetwork{typeof(asys),E,V}(asys, reactions, cuts, halves,
                                            Vector{Float64}(undef, n), Vector{Float64}(undef, n),
                                            Vector{Bool}(undef, n))
    return rate!(net; kwargs...)
end

"""
    rate!(net::ReactionNetwork; fwdkernel=Returns(1.0), bwdkernel=fwdkernel)

Set the rates of the reactions of `net` in place. The kernels should take a [`Reaction`](@ref)
and return a (non-negative) rate.
"""
function rate!(net::ReactionNetwork; fwdkernel=Returns(1.0), bwdkernel=fwdkernel)
    for r in eachindex(net._reactions)
        rxn = Reaction(net, r)
        net._fwdrates[r] = fwdkernel(rxn)
        net._bwdrates[r] = bwdkernel(rxn)
        net.active[r] = net._fwdrates[r] != 0 || net._bwdrates[r] != 0
    end
    return net
end

Base.length(net::ReactionNetwork) = length(net._reactions)
Base.eltype(::Type{N}) where {N<:ReactionNetwork} = Reaction{N}
Base.getindex(net::ReactionNetwork, r::Integer) = Reaction(net, r)
Base.iterate(net::ReactionNetwork, r=1) = r > length(net) ? nothing : (net[r], r + 1)

nstructures(net::ReactionNetwork) = nstructures(net.assemblysystem)
nreactions(net::ReactionNetwork) = length(net)

function Base.show(io::Core.IO, net::ReactionNetwork)
    return print(io, "ReactionNetwork[$(count(net.active))/$(length(net)) reactions active]")
end
