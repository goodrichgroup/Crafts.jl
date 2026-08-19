module CraftsJuMPExt

using Crafts, JuMP, SparseArrays, LinearAlgebra

# Snap a float direction to a small exact rational; the true generators of count cones have small
# denominators, so this usually recovers them exactly (but nothing certifies that).
_snap(v) = rationalize.(BigInt, v ./ maximum(abs, v); tol=1e-8)

# `A[rows, :] == -I` without materializing a comparison: exactly one entry per row/column, every
# diagonal entry -1.
function _isnegidentity(N)
    size(N, 1) == size(N, 2) || return false
    n = size(N, 1)
    count(!iszero, N) == n || return false
    return all(N[i, i] == -1 for i in 1:n)
end

function Crafts._lporacle(optimizer::Union{Type,Function,MOI.OptimizerWithAttributes},
                          P, A, linearity)
    m, n = size(A)
    islin = falses(m)
    islin[linearity] .= true
    N = A[.!islin, :]
    _isnegidentity(N) && return _cgoracle(optimizer, P, A[islin, :])

    # Generic fallback: the whole system as one model, built once; each probe only swaps the
    # objective, so warm starts carry over between solves. The probed directions must be exact
    # facet candidates: over a cone the support value is 0 or +∞, discontinuous in the direction,
    # so a direction tilted off a facet by float noise would read as unbounded no matter the
    # tolerance.
    Af = sparse(Float64.(A))
    Pf = Float64.(P)
    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[1:n])
    @constraint(model, Af[.!islin, :] * x .<= 0)
    @constraint(model, Af[islin, :] * x .== 0)

    return function (c)
        # normalize exactly before the float conversion: raw facet normals can carry huge integer
        # entries whose Float64 images overflow or drown the solver tolerances
        w = Float64.(transpose(P) * (c ./ maximum(abs, c)))
        @objective(model, Max, w' * x)
        optimize!(model)
        st = termination_status(model)
        if st == MOI.OPTIMAL
            return true, Vector{Rational{BigInt}}[]
        elseif primal_status(model) == MOI.INFEASIBILITY_CERTIFICATE
            # unbounded: the certificate is an improving ray of the source cone
            return false, [_snap(Pf * value.(x))]
        end
        error("the LP solver returned `$st` on a support probe; the cone cannot be trusted")
    end
end

# Cached-column column generation for the structured cone `{μ >= 0 : L*μ = 0}`. Every support LP
# has a basic solution touching at most `size(L, 1) + 1` of the (possibly hundreds of thousands
# of) columns, so the restricted master only ever holds a small working set, shared across
# probes; the remaining columns exist only as a cached sparse matrix swept by the pricing step.
# The master proves `sup = 0` on the working set, the pricing sweep extends the proof to every
# column, and an unbounded master hands back the improving ray.
function _cgoracle(optimizer, P, L)
    nr, n = size(L)
    Lf = sparse(Float64.(L))
    LfT = sparse(transpose(Lf))
    PfT = permutedims(Float64.(P))

    model = Model(optimizer)
    set_silent(model)
    # the master is small; without presolve the solver can always hand back the unboundedness
    # certificate (with it, HiGHS sometimes fails the solve outright instead)
    try
        set_attribute(model, "presolve", "off")
    catch
    end
    bal = @constraint(model, [r = 1:nr], zero(AffExpr) == 0)

    S = Int[]                 # columns in the master, in variable order
    vars = VariableRef[]
    added = falses(n)
    rows = rowvals(Lf)
    vals = nonzeros(Lf)
    function addcolumn!(j)
        added[j] && return false
        v = @variable(model, lower_bound = 0)
        push!(S, j)
        push!(vars, v)
        added[j] = true
        for k in nzrange(Lf, j)
            set_normalized_coefficient(bal[rows[k]], v, vals[k])
        end
        return true
    end

    # JuMP's dual sign convention for equality rows is decided per solve by dual feasibility on
    # the working set: under the correct sign, no included column has positive reduced cost.
    function _lambda(w)
        λ = dual.(bal)
        excess(s) = maximum(w[j] - s * dot(view(Lf, :, j), λ) for j in S)
        return excess(1) <= excess(-1) ? λ : -λ
    end

    batch = 100
    return function (c)
        # exact normalization first — see the generic oracle
        w = PfT * Float64.(c ./ maximum(abs, c))
        tol = 1e-7 * (1 + maximum(abs, w))
        for _ in 1:100_000
            λ = zeros(nr)
            if !isempty(vars)
                @objective(model, Max, sum(w[S[k]] * vars[k] for k in eachindex(vars)))
                optimize!(model)
                st = termination_status(model)
                if primal_status(model) == MOI.INFEASIBILITY_CERTIFICATE
                    ray = value.(vars)
                    point = zeros(size(PfT, 2))
                    for (k, j) in enumerate(S)
                        iszero(ray[k]) || (point .+= ray[k] .* view(PfT, j, :))
                    end
                    return false, [_snap(point)]
                end
                st == MOI.OPTIMAL ||
                    error("the LP solver returned `$st` (\"$(raw_status(model))\") on a support " *
                          "probe with $(length(vars)) master columns; the cone cannot be trusted")
                λ = _lambda(w)
            end

            # pricing: reduced costs of every cached column against the master's duals
            red = w .- LfT * λ
            order = partialsortperm(red, 1:min(batch, n); rev=true)
            grew = false
            for j in order
                red[j] > tol || break
                grew |= addcolumn!(j)
            end
            grew || return true, Vector{Rational{BigInt}}[]
        end
        error("column generation did not converge; the support probe cannot be trusted")
    end
end

# Float sweep of the wedge slab `{μ ≥ 0 : Lμ = 0, (hᵀP)μ = 1}`: which environments can carry
# weight in a count vector beyond the survivor-cone facet `h`. Over-collecting is sound (any
# superset of the true support keeps the wedge restriction lossless), so borderline columns are
# included: unbounded or positive-within-tolerance both count. `h` is normalized exactly before
# the float conversion — survivor-hull facet normals carry huge integer entries.
function Crafts._wedgemembers(optimizer::Union{Type,Function,MOI.OptimizerWithAttributes},
                              L, P, h)
    n = size(L, 2)
    hn = h ./ maximum(abs, h)
    w = Float64.(vec(transpose(Rational{BigInt}.(hn)) * P))
    Lf = sparse(Float64.(L))

    model = Model(optimizer)
    set_silent(model)
    @variable(model, μ[1:n] >= 0)
    @constraint(model, Lf * μ .== 0)
    @constraint(model, w' * μ == 1)

    members = falses(n)
    undecided = trues(n)
    for j in 1:n
        undecided[j] || continue
        @objective(model, Max, μ[j])
        optimize!(model)
        st = termination_status(model)
        if st == MOI.OPTIMAL
            x = value.(μ)
            for k in 1:n
                x[k] > 1e-9 || continue
                members[k] = true
                undecided[k] = false
            end
            x[j] > 1e-9 || (undecided[j] = false)
        elseif st == MOI.INFEASIBLE
            return Int[]   # empty slab: nothing lies beyond this facet
        else
            # unbounded or ambiguous: include conservatively
            members[j] = true
            undecided[j] = false
        end
    end
    return findall(members)
end

# Support of the fiber `{μ ≥ 0 : Sμ ≤ 0, Lμ = 0, Pμ = m}`: every environment that can carry weight
# in a count vector realizing `m`. One model, one objective swap per environment — the exact route
# spends an lrs subprocess per environment instead, which is what dominates certification.
#
# Over-collecting is sound (the support only decides which environments get refined), so an
# ambiguous solve includes the environment rather than certifying anything. `m` is normalized: the
# fiber scales with it, so the support is unchanged and the right-hand side stays O(1).
function Crafts._fibermembers(optimizer::Union{Type,Function,MOI.OptimizerWithAttributes},
                              S, L, P, m)
    n = size(L, 2)
    mr = Rational{BigInt}.(collect(m))
    mf = Float64.(mr ./ maximum(abs, mr))
    model = Model(optimizer)
    set_silent(model)
    @variable(model, μ[1:n] >= 0)
    @constraint(model, sparse(Float64.(L)) * μ .== 0)
    size(S, 1) == 0 || @constraint(model, sparse(Float64.(S)) * μ .<= 0)
    @constraint(model, sparse(Float64.(P)) * μ .== mf)

    members = falses(n)
    undecided = trues(n)
    for j in 1:n
        undecided[j] || continue
        @objective(model, Max, μ[j])
        optimize!(model)
        st = termination_status(model)
        if st == MOI.OPTIMAL
            x = value.(μ)
            for k in 1:n
                x[k] > 1e-9 || continue
                members[k] = true
                undecided[k] = false
            end
            x[j] > 1e-9 || (undecided[j] = false)
        elseif st == MOI.INFEASIBLE
            return Int[]   # empty fiber: `m` is not in the cone at all
        else
            members[j] = true
            undecided[j] = false
        end
    end
    return findall(members)
end

# Float feasibility of the fiber `{μ ≥ 0 : Lμ = 0, Pμ = m}`, with the one verdict that matters
# — infeasibility, i.e. exact elimination of the direction `m` — accepted only after its Farkas
# certificate `Lᵀλ + Pᵀν ≤ 0, νᵀm > 0` verifies in exact rational arithmetic. Feasibility is
# reported as a snapped count vector without certification (the safe direction), and an
# unverifiable certificate comes back as `missing`.
function Crafts._fiberfeasible(optimizer::Union{Type,Function,MOI.OptimizerWithAttributes},
                               S, L, P, m)
    n = size(L, 2)
    Sf = sparse(Float64.(S))
    Lf = sparse(Float64.(L))
    Pf = Float64.(P)

    model = Model(optimizer)
    set_silent(model)
    try
        set_attribute(model, "presolve", "off")   # keeps the dual ray available
    catch
    end
    @variable(model, μ[1:n] >= 0)
    slk = @constraint(model, Sf * μ .<= 0)
    bal = @constraint(model, Lf * μ .== 0)
    proj = @constraint(model, Pf * μ .== Float64.(m))
    optimize!(model)

    if termination_status(model) == MOI.OPTIMAL
        return rationalize.(BigInt, value.(μ); tol=1e-8)
    end
    dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE || return missing
    σ = rationalize.(BigInt, dual.(slk); tol=1e-10)
    λ = rationalize.(BigInt, dual.(bal); tol=1e-10)
    ν = rationalize.(BigInt, dual.(proj); tol=1e-10)
    y = transpose(Rational{BigInt}.(S)) * σ + transpose(Rational{BigInt}.(L)) * λ +
        transpose(Rational{BigInt}.(P)) * ν
    s = transpose(ν) * Rational{BigInt}.(m)
    # either sign orientation of the ray proves emptiness; the slack duals must carry the sign
    # matching their inequality side
    ((all(<=(0), y) && s > 0 && all(>=(0), σ)) ||
     (all(>=(0), y) && s < 0 && all(<=(0), σ))) && return nothing
    return missing
end

end
