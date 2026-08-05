"""
    simulate_kinetics(asys::AssemblySystem, ξ; T, kernel, [initial_densities=nothing, saveat=[]])

Simulate the assembly kinetics of the structures of `asys` as a function of
`ξ`, a vector containing chemical potentials and binding energies. Unless `initial_densities` is specified,
the initial state consists of monomers at the concentrations implied by `ξ`.

Keyword arguments:
- `T`: total simulation time.
- `kernel`: aggregation kernel used to build the reaction network.
- `initial_densities` (optional): initial densities of all structures at t=0.
- `saveat` (optional): time points at which the simulation state should be stored.

Returns a vector of time points and a matrix of structure number densities of shape `(nstructures, ntimepoints)`.
"""
function simulate_kinetics(asys::AssemblySystem, ξ; T, initial_densities=nothing, saveat=[], kernel, kwargs...)
    _check_parameterlength(asys, ξ)
    ts, ρs = _simulate_kinetics(asys, ξ; Ts=(0,T), kernel, initial_densities, saveat, kwargs...)
    return ts, ρs
end
function simulate_kinetics(asys::AssemblySystem, ϕs, εs; kwargs...)
    _check_parameterlength(asys, ϕs, εs)
    return simulate_kinetics(asys, [chemicalpotentials(asys, ϕs, εs); εs]; kwargs...)
end


function cleave(anatomy::AbstractGraph, es::Vector{<:AbstractEdge})
    anatomy = copy(anatomy)

    for edge in es
        revedge = reverse(edge)
        if revedge ∉ edges(anatomy)
            error("Cannot cleave a nondirected (interior) edge.")
        end
        rem_edge!(anatomy, edge)
        rem_edge!(anatomy, revedge)
    end

    comps = connected_components(anatomy) 
    return anatomy, comps
end
cleave(anatomy::AbstractGraph, edge::AbstractEdge) = cleave(anatomy, [edge])

function listbonds(bondcounts)
    # Convert bond counts, i.e. [0, 1, 0, 0, 2]
    # into a list of bonds (with repeats), i.e. [2, 5, 5]
    bondlist = zeros(Int, sum(bondcounts))

    k = 1
    for (i, b) in pairs(bondcounts)
        if b > 0
            bondlist[k:k+b-1] .= i
            k += b
        end
    end
    return bondlist
end

function generate_reactions(g; maxlevel)
    cuts, components = generate_cuts(g; maxlevel)
    # Reaction stores graph, graph -> graph at this point
    reactions = [(g[c1], g[c2], g) for (c1, c2) in components]
    return reactions, cuts, components
end


function are_separated(vi, vj, vs1, vs2)
    return (vi ∈ vs1 && vj ∈ vs2) || 
           (vi ∈ vs2 && vj ∈ vs1)
end

function generate_cuts(g; maxlevel)
    edges = collect(Roly.exterior_edges(g))
    cuts = Vector{eltype(edges)}[]
    halfs = Vector{Vector{Int32}}[]
    isempty(edges) && return cuts, halfs

    current_cut = [first(edges)]

    # Perform a depth-first backtracking search to 
    # generate all possible cuts of length <= maxlevel
    while true
        current_edge = current_cut[end]

        # Remove all edges of the current cut from the input graph and return connected components
        _, component_vertices = cleave(g, current_cut)

        if length(component_vertices) > 1
            # Cut was successful in separating the graph

            # Check if all vertices of the cut are in different components
            # If that is the case, add the cut to the output
            srcs = [e.src for e in current_cut]
            dsts = [e.dst for e in current_cut]
            if all(x -> are_separated(x..., component_vertices...), zip(srcs, dsts))
                push!(cuts, copy(current_cut))
                push!(halfs, component_vertices)
            end

            nextedge_idx = nothing # initiate an upward traverse
        elseif length(current_cut) == maxlevel
            nextedge_idx = nothing # initiate an upward traverse
        else
            # Cut is incomplete, keep traversing downward
            nextedge_idx = findfirst(>(current_edge), edges) # TODO use searchsorted
        end

        # If no additional edge can be added to the cut, traverse upward until a viable
        # branch is found
        while isnothing(nextedge_idx)
            isempty(current_cut) && @goto finished
            current_edge = pop!(current_cut)
            nextedge_idx = findfirst(>(current_edge), edges) # TODO use searchsorted
        end

        next_edge = edges[nextedge_idx]
        push!(current_cut, next_edge)
    end

    @label finished

    return cuts, halfs
end

function make_kernels(reactions, cuts, components, aggkernel=nothing, brkkernel=nothing)
    if isnothing(aggkernel)
        aggkernel = (k, cut, component) -> 1
    end
    if isnothing(brkkernel)
        brkkernel = (k, cut, component) -> 1
    end

    ks = [aggkernel(k, cut, component) for ((_, _, k), cut, component) in zip(reactions, cuts, components)]
    fs = [brkkernel(k, cut, component) for ((_, _, k), cut, component) in zip(reactions, cuts, components)]
    return ks, fs
end

function generate_reactionnetwork(strs; maxlevel, aggkernel=nothing, brkkernel=aggkernel)
    gs = [s.graphrep for s in strs]
    foreach(canonize!, gs)
    ids = Dict(hash(g)=>i for (i, g) in enumerate(gs))

    reactions = NTuple{3,Int}[] # reactions in the form i, j <--> k
    cuts = []
    components = []

    for g in gs
        _cuts, _components = generate_cuts(g; maxlevel)
        graph_reacts = [(g[c1], g[c2], g) for (c1, c2) in _components]
        foreach(gr->map(canonize!, gr), graph_reacts)

        isempty(graph_reacts) && continue
        
        for (cut, component, greact) in zip(_cuts, _components, graph_reacts)
            react = map(x->ids[hash(x)], greact)

            push!(reactions, react)
            push!(cuts, cut)
            push!(components, component)
        end
    end
    ks, fs = make_kernels(reactions, cuts, components, aggkernel, brkkernel)

    nonzero_rs = filter(r->ks[r] ≉ 0 || fs[r] ≉ 0, eachindex(reactions))
    reactions = reactions[nonzero_rs]
    ks = ks[nonzero_rs]
    fs = fs[nonzero_rs]

    return reactions, ks, fs
end

function kinetic_network(asys::AssemblySystem, ξ; maxbonds, kernel, brkkernel=kernel)
    reactions, ks, fs = generate_reactionnetwork(structures(asys); maxlevel=maxbonds, aggkernel=kernel, brkkernel)

    log_ρs = logdensities(asys, ξ)
    fs = [fs[r] * exp(log_ρs[i] + log_ρs[j] - log_ρs[k]) for (r, (i,j,k)) in enumerate(reactions)]
    kscale = exp(mean(log.(fs))) # geometrical mean
    fs /= kscale

    function update_step!(du, u, p, t)
        du .= 0

        for r in eachindex(reactions)
            i, j, k = reactions[r]

            Ka = ks[r] * u[i] * u[j]
            Kb = fs[r] * u[k]

            Rij = Kb - Ka # We don't need a factor of 2 for i == j, because we add it twice in that case!

            du[i] += Rij
            du[j] += Rij
            du[k] -= Rij
        end
        return
    end
    return update_step!, kscale
end

function _simulate_kinetics(asys::AssemblySystem, ξ; Ts, kernel, brkkernel=kernel, maxbonds=Inf, initial_densities=nothing, alg=Rodas5(), saveat=[], abstol=nothing, solver_kwargs...)
    np = nspecies(asys)
    nstr = nstructures(asys)

    step, kscale = kinetic_network(asys, ξ; kernel, brkkernel, maxbonds)
    tscale = inv(kscale)
    ρscale = kscale

    if isnothing(initial_densities)
        initial_densities = vcat(particledensities(asys, ξ), zeros(nstr - np))
    end
    if isnothing(abstol)
        abstol = sum(initial_densities) * 1e-8
    end
    
    initial_densities /= ρscale
    Ts = Ts ./ tscale
    saveat = saveat ./ tscale
    abstol /= ρscale

    prob = ODEProblem(step, initial_densities, Ts)
    sol = solve(prob, alg; saveat, abstol, solver_kwargs...)

    ts = sol.t * tscale
    us = reduce(hcat, sol.u * ρscale)
    return ts, us
end

function stability_matrix(asys::AssemblySystem; symmetrize=false, kernel, brkkernel=kernel, maxbonds)
    reactions, ks, fs = generate_reactionnetwork(structures(asys); maxlevel=maxbonds, aggkernel=kernel, brkkernel)
    M = compositionmatrix(asys)
    Ωs = partitionfunctions(asys)

    function Sfn!(S, ξ)
        S .= 0
        ρeq = densities(ξ, M, Ωs)

        for r in eachindex(reactions)
            i, j, k = reactions[r]

            Ka = ks[r]
            Kb = (ρeq[i] * ρeq[j] / ρeq[k]) * fs[r]

            S[i, i] += -Ka * ρeq[j] 
            S[i, j] += -Ka * ρeq[i] 
            S[i, k] += Kb

            S[j, i] += -Ka * ρeq[j] 
            S[j, j] += -Ka * ρeq[i] 
            S[j, k] += Kb

            S[k, i] += Ka * ρeq[j] 
            S[k, j] += Ka * ρeq[i] 
            S[k, k] += -Kb
        end
        return S 
    end

    function ∂S∂μfn!(S, ξ)
        S .= 0
        ρeq = densities(ξ, M, Ωs)
        ∂ρeq = permutedims(∂ρ∂μ(ξ, M, Ωs))

        scratch = zero(ξ)

        for r in eachindex(reactions)
            i, j, k = reactions[r]

            Ka = ks[r]
            Kb = fs[r] 

            @views @. begin 
                scratch = -Ka * ∂ρeq[:, j] 
                S[i, i, :] += scratch
                S[j, i, :] += scratch
                S[k, i, :] -= scratch


                scratch = -Ka * ∂ρeq[:, i] 
                S[i, j, :] += scratch
                S[j, j, :] += scratch
                S[k, j, :] -= scratch

                scratch = Kb *(-ρeq[i] * ρeq[j] / ρeq[k]^2 * ∂ρeq[:, k] + 
                            ∂ρeq[:, i] * ρeq[j] / ρeq[k] 
                            + ∂ρeq[:, j] * ρeq[i] / ρeq[k])
                S[i, k, :] += scratch
                S[j, k, :] += scratch
                S[k, k, :] -= scratch
            end
        end
        return S
    end

    #####################
    #### Return inv(D) S D, where D = sqrt(ρᵢ)
    function Ssym_fn!(S, ξ; scale=1)
        ρeq = densities(ξ, M, Ωs)
        ρeq_sqrt = sqrt.(ρeq)

        S .= 0
        for r in eachindex(reactions)
            i, j, k = reactions[r]

            Ka = ks[r] * scale
            Kb = fs[r] * scale

            ij = ρeq_sqrt[i] * ρeq_sqrt[j] 
            ik = ρeq_sqrt[i] / ρeq_sqrt[k] * ρeq[j]
            jk = ρeq_sqrt[j] / ρeq_sqrt[k] * ρeq[i]

            S[i, i] += -Ka * ρeq[j]
            S[i, j] += -Ka * ij
            S[i, k] += Kb * ik

            S[j, i] += -Ka * ij
            S[j, j] += -Ka * ρeq[i]
            S[j, k] += Kb * jk

            S[k, i] += Ka * ik
            S[k, j] += Ka * jk
            S[k, k] += -Kb * (ρeq[i] * ρeq[j] / ρeq[k])
        end
        return S
    end

    function ∂Ssym∂μ_fn!(S, ξ; scale=1)
        ρeq = densities(ξ, M, Ωs)
        ρeq_sqrt = sqrt.(ρeq)
        ∂ρeq = permutedims(∂ρ∂μ(ξ, M, Ωs))

        S .= 0
        scratch = zero(ξ)

        for r in eachindex(reactions)
            i, j, k = reactions[r]

            Ka = ks[r] * scale
            Kb = fs[r] * scale

            @views @. begin
                scratch = -Ka * ((ρeq_sqrt[j] / ρeq_sqrt[i]) * ∂ρeq[:, i] + (ρeq_sqrt[i] / ρeq_sqrt[j]) * ∂ρeq[:, j]) / 2
                S[i, j, :] += scratch
                S[j, i, :] += scratch

                # ik
                scratch = (ρeq_sqrt[i] / ρeq_sqrt[k] * ∂ρeq[:, j] + ρeq[j] / (2ρeq_sqrt[i] * ρeq_sqrt[k]) * ∂ρeq[:, i] -
                    ρeq_sqrt[i] * ρeq[j] / (2ρeq_sqrt[k]^3) * ∂ρeq[:, k])
                S[i, k, :] += Kb * scratch
                S[k, i, :] += Ka * scratch

                # jk
                scratch = (ρeq_sqrt[j] / ρeq_sqrt[k] * ∂ρeq[:, i] + ρeq[i] / (2ρeq_sqrt[j] * ρeq_sqrt[k]) * ∂ρeq[:, j] -
                ρeq_sqrt[j] * ρeq[i] / (2ρeq_sqrt[k]^3) * ∂ρeq[:, k])

                S[j, k, :] += Kb * scratch
                S[k, j, :] += Ka * scratch

                S[i, i, :] += -Ka * ∂ρeq[:, j] 
                S[j, j, :] += -Ka * ∂ρeq[:, i] 
                S[k, k, :] += -Kb * (-ρeq[i] * ρeq[j] / ρeq[k]^2 * ∂ρeq[:, k] + ∂ρeq[:, i] * ρeq[j] / ρeq[k] + ∂ρeq[:, j] * ρeq[i] / ρeq[k])
            end
        end
        return S
    end

    if symmetrize
        return Ssym_fn!, ∂Ssym∂μ_fn!
    else
        return Sfn!, ∂S∂μfn!
    end
end

function τc(S; np)
    ns = size(S, 1)

    C = maximum(abs, S)
    S = S / C

    λ = partialsort!(schur(S).values * C, ns-np, by=real)
    return -inv(real(λ))
end
