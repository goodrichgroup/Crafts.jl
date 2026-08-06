"""
    AssemblySystem

Wraps a set of Roly.jl `BindingRules` together with an [`EntropyModel`](@ref) used to compute structure
partition functions. Structure enumeration, the composition matrix, and partition functions are
computed lazily and cached on first use.
"""
mutable struct AssemblySystem{BR<:BindingRules,EM<:EntropyModel,P<:Polyform}
    const bindingrules::BR
    const entropymodel::EM
    const maxsize::Float64
    const verbose::Bool

    structures::Union{Nothing,Vector{P}}
    M::Union{Nothing,Matrix{Int}}
    Ωs::Union{Nothing,Vector{Float64}}
    status::Union{Nothing,RSStatus}
end

"""
    AssemblySystem(bindingrules::BindingRules, entropymodel::EntropyModel; maxsize=Inf, verbose=true)

Construct an `AssemblySystem` from `bindingrules` and an [`EntropyModel`](@ref) used to compute structure
partition functions. Only structures of at most `maxsize` particles are considered.
"""
function AssemblySystem(bindingrules::BindingRules, entropymodel::EntropyModel; maxsize=Inf,
                        verbose=true)
    (maxsize > 0 && (isinf(maxsize) || isinteger(maxsize))) ||
        throw(ArgumentError("maxsize=$maxsize must be a positive integer or Inf."))
    checkpotential(entropymodel.potential, bindingrules,
                   entropymodel.embed3d ? 3 : dimension(bindingrules))

    P = typeof(Polyform(bindingrules))
    return AssemblySystem{typeof(bindingrules),typeof(entropymodel),P}(
        bindingrules, entropymodel, maxsize, verbose, nothing, nothing, nothing, nothing,
    )
end

function Base.show(io::Core.IO, asys::AssemblySystem)
    cap = isfinite(asys.maxsize) ? ", maxsize=$(Int(asys.maxsize))" : ""
    cut = asys.status == MaxDepthReached ? ", truncated" : ""
    return print(io, "AssemblySystem[n=$(nspecies(asys)), k=$(nbonds(asys))$cap$cut]")
end

nspecies(asys::AssemblySystem) = nspecies(asys.bindingrules)
nbonds(asys::AssemblySystem) = nbonds(asys.bindingrules)
dimension(asys::AssemblySystem) = dimension(asys.bindingrules)

"""
    iscomplete(asys::AssemblySystem)

Whether `asys` covers every structure its rules allow, i.e. whether `maxsize` cut the enumeration short.
`missing` if nothing has been enumerated yet.
"""
iscomplete(asys::AssemblySystem) = asys.status === nothing ? missing : asys.status == Finished

_iscached(asys::AssemblySystem, kw) = kw.maxsize == asys.maxsize && length(kw) == 1

# raise a warning if the enumeration is truncated a a given size
function _warnifincomplete(asys::AssemblySystem)
    (asys.verbose && asys.status != Finished) || return nothing
    cap = isfinite(asys.maxsize) ? Int(asys.maxsize) : asys.maxsize
    @warn "AssemblySystem enumeration stopped at $(asys.status) with maxsize=$cap, so these are " *
          "not all the structures the rules allow. Everything computed from them is a truncation; " *
          "see `iscomplete`."
    return nothing
end

function _tomatrix(asys::AssemblySystem, rows)
    M = Matrix{Int}(undef, length(rows), nspecies(asys) + nbonds(asys))
    for (i, row) in enumerate(rows)
        M[i, :] .= row
    end
    return M
end

function _enumerate_structures(asys::AssemblySystem{BR,F,P}, kwargs) where {BR,F,P}
    strs, sizes = P[], Int[]
    res = polyenum(asys.bindingrules; kwargs...) do s, n
        push!(strs, copy(s))
        push!(sizes, n)
        ACCEPT
    end
    return strs[sortperm(sizes; alg=MergeSort)], res.status
end

function _enumerate_compositions(asys::AssemblySystem, kwargs)
    rows, sizes = Vector{Int}[], Int[]
    res = polyenum(asys.bindingrules; kwargs...) do s, n
        push!(rows, composition(s))
        push!(sizes, n)
        ACCEPT
    end
    return _tomatrix(asys, (rows[i] for i in sortperm(sizes; alg=MergeSort))), res.status
end

_entropies(asys::AssemblySystem, strs) = Float64[entropy(asys.entropymodel, p) for p in strs]

"""
    structures(asys::AssemblySystem; enum_kwargs...)

Return the enumerated structures of `asys`, computing and caching them on first use.
"""
function structures(asys::AssemblySystem; maxsize=asys.maxsize, kwargs...)
    kw = (; maxsize, kwargs...)
    _iscached(asys, kw) || return first(_enumerate_structures(asys, kw))

    strs = asys.structures
    if strs === nothing
        firstlook = asys.status === nothing
        strs, asys.status = _enumerate_structures(asys, kw)
        firstlook && _warnifincomplete(asys)
        asys.structures = strs
        # Derive `M` from this exact ordering, so it can never disagree with `Ωs`.
        asys.M = _tomatrix(asys, (composition(s) for s in strs))
    end
    return strs
end

"""
    compositionmatrix(asys::AssemblySystem; enum_kwargs...)

Return the composition matrix of `asys`, computing and caching it on first use. Unless the structures
themselves are already cached, they are not retained.
"""
function compositionmatrix(asys::AssemblySystem; maxsize=asys.maxsize, kwargs...)
    kw = (; maxsize, kwargs...)
    _iscached(asys, kw) || return first(_enumerate_compositions(asys, kw))

    M = asys.M
    if M === nothing
        firstlook = asys.status === nothing
        M, asys.status = _enumerate_compositions(asys, kw)
        firstlook && _warnifincomplete(asys)
        asys.M = M
    end
    return M
end

"""
    nstructures(asys::AssemblySystem; enum_kwargs...)

Return the number of enumerated structures of `asys`, without retaining the structures themselves.
"""
nstructures(asys::AssemblySystem; kwargs...) = size(compositionmatrix(asys; kwargs...), 1)

"""
    countstructures(asys::AssemblySystem; kwargs...)

Count how many structures `asys` *would* enumerate, without enumerating them. Cheap enough to run
before committing to [`structures`](@ref) or [`partitionfunctions`](@ref) on a large system.

Returns a `StructureCount`, exact when the enumeration fits in `exact_budget` and a statistical estimate
otherwise — as opposed to [`nstructures`](@ref), which returns the plain number of structures actually
enumerated and pays for the enumeration to do it.
"""
countstructures(asys::AssemblySystem; kwargs...) = Roly.countpolyforms(asys.bindingrules; kwargs...)

"""
    partitionfunctions(asys::AssemblySystem; enum_kwargs...)

Return the per-structure partition functions of `asys`, computed via [`entropy`](@ref) using
`asys.bondpotential`, computing and caching them on first use.
"""
function partitionfunctions(asys::AssemblySystem; maxsize=asys.maxsize, kwargs...)
    _iscached(asys, (; maxsize, kwargs...)) || return _entropies(asys, structures(asys; maxsize, kwargs...))

    Ωs = asys.Ωs
    if Ωs === nothing
        Ωs = _entropies(asys, structures(asys))
        asys.Ωs = Ωs
    end
    return Ωs
end
