module CraftsJuMPExt

using Crafts, JuMP, SparseArrays, LinearAlgebra

# Snap a float direction to a small exact rational; the true generators of count cones have small
# denominators, so this usually recovers them exactly (but nothing certifies that).
_snap(v) = rationalize.(BigInt, v ./ maximum(abs, v); tol=1e-8)

# `N == -I` without materializing a comparison: exactly one entry per row/column, every diagonal
# entry -1.
function _isnegidentity(N)
    size(N, 1) == size(N, 2) || return false
    n = size(N, 1)
    count(!iszero, N) == n || return false
    return all(N[i, i] == -1 for i in 1:n)
end

function Crafts._lporacle(optimizer::Union{Type,Function,MOI.OptimizerWithAttributes},
                          P, A, linearity; verbose::Bool=false)
    m, n = size(A)
    islin = falses(m)
    islin[linearity] .= true
    N = A[.!islin, :]
    _isnegidentity(N) && return _ipmoracle(optimizer, P, A[islin, :]; verbose)

    # Generic fallback for an arbitrary cone: the whole system as one model, built once; each
    # probe only swaps the objective. The probed directions must be exact facet candidates: over
    # a cone the support value is 0 or +∞, discontinuous in the direction, so a direction tilted
    # off a facet by float noise would read as unbounded no matter the tolerance.
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

# The structured cone `{μ >= 0 : L*μ = 0}` as one sparse model, sliced by `Σμ = 1` so that every
# support LP is bounded: its optimum is 0 iff the probed facet is valid, and otherwise the
# maximizer itself is an improving point. Built once, the objective swapped per probe, solved by
# an interior point method with crossover (HiGHS options; silently ignored by other solvers).
# Simplex, warm-started or not, grinds on the degenerate optimal faces of these LPs (a valid
# facet can be tight on half of all environments); interior point is indifferent to their
# dimension.
#
# A valid verdict is checked exactly where the duals allow: `λ` from the balance rows,
# rationalized, certifies the facet iff `Pᵀc - Lᵀλ <= 0` holds in exact arithmetic. Verdicts that
# fail the check are still accepted (there is no exact fallback at this size) and counted.
function _ipmoracle(optimizer, P, L; verbose::Bool=false)
    n = size(L, 2)
    Lf = sparse(Float64.(L))
    LT = transpose(Rational{BigInt}.(L))
    PT = transpose(Rational{BigInt}.(P))
    PfT = permutedims(Float64.(P))

    model = Model(optimizer)
    set_silent(model)
    for (k, v) in (("solver", "ipm"), ("run_crossover", "on"), ("presolve", "on"))
        try
            set_attribute(model, k, v)
        catch
        end
    end
    @variable(model, μ[1:n] >= 0)
    bal = @constraint(model, Lf * μ .== 0)
    @constraint(model, sum(μ) == 1)

    uncertified = 0
    return function (c)
        cn = c ./ maximum(abs, c)
        w = PfT * Float64.(cn)
        tol = 1e-7 * (1 + maximum(abs, w))
        @objective(model, Max, w' * μ)
        t = @elapsed optimize!(model)
        st = termination_status(model)
        st == MOI.OPTIMAL ||
            error("the LP solver returned `$st` (\"$(raw_status(model))\") on a support probe; " *
                  "the cone cannot be trusted")
        x = value.(μ)
        val = objective_value(model)
        if val > tol
            verbose && (println("    ipm: violated by ", round(val; sigdigits=3), ", support ",
                                count(>(1e-9), x), " (", round(t; digits=1), "s)"); flush(stdout))
            return false, [_snap(PfT' * x)]
        end
        λ = rationalize.(BigInt, dual.(bal); tol=1e-10)
        exact = all(<=(0), PT * Rational{BigInt}.(cn) - LT * λ)
        exact || (uncertified += 1)
        verbose && (println("    ipm: valid, ", exact ? "certified exactly" : "float only (#$uncertified)",
                            " (", round(t; digits=1), "s)"); flush(stdout))
        return true, Vector{Rational{BigInt}}[]
    end
end

# Support of the fiber `{μ ≥ 0 : Lμ = 0, Pμ = m}`: every environment that can carry weight in a
# count vector realizing `m`. One model, one objective swap per environment; the exact route
# spends an exact LP per environment instead, which is what dominates certification.
#
# Over-collecting is sound (the support only decides which environments get refined), so an
# ambiguous solve includes the environment rather than certifying anything. `m` is normalized: the
# fiber scales with it, so the support is unchanged and the right-hand side stays O(1).
function Crafts._fibermembers(optimizer::Union{Type,Function,MOI.OptimizerWithAttributes},
                              L, P, m)
    n = size(L, 2)
    mr = Rational{BigInt}.(collect(m))
    mf = Float64.(mr ./ maximum(abs, mr))
    model = Model(optimizer)
    set_silent(model)
    @variable(model, μ[1:n] >= 0)
    @constraint(model, sparse(Float64.(L)) * μ .== 0)
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

# Float feasibility of the fiber `{μ ≥ 0 : Lμ = 0, Pμ = m}`, with the one verdict that matters,
# infeasibility, i.e. exact elimination of the direction `m`, accepted only after its Farkas
# certificate `Lᵀλ + Pᵀν ≤ 0, νᵀm > 0` verifies in exact rational arithmetic. Feasibility is
# reported as a snapped count vector without certification (the safe direction), and an
# unverifiable certificate comes back as `missing`.
function Crafts._fiberfeasible(optimizer::Union{Type,Function,MOI.OptimizerWithAttributes},
                               L, P, m)
    n = size(L, 2)
    Lf = sparse(Float64.(L))
    Pf = Float64.(P)

    model = Model(optimizer)
    set_silent(model)
    try
        set_attribute(model, "presolve", "off")   # keeps the dual ray available
    catch
    end
    @variable(model, μ[1:n] >= 0)
    bal = @constraint(model, Lf * μ .== 0)
    proj = @constraint(model, Pf * μ .== Float64.(m))
    optimize!(model)

    if termination_status(model) == MOI.OPTIMAL
        return rationalize.(BigInt, value.(μ); tol=1e-8)
    end
    dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE || return missing
    λ = rationalize.(BigInt, dual.(bal); tol=1e-10)
    ν = rationalize.(BigInt, dual.(proj); tol=1e-10)
    y = transpose(Rational{BigInt}.(L)) * λ + transpose(Rational{BigInt}.(P)) * ν
    s = transpose(ν) * Rational{BigInt}.(m)
    # either sign orientation of the ray proves emptiness
    ((all(<=(0), y) && s > 0) || (all(>=(0), y) && s < 0)) && return nothing
    return missing
end

end
