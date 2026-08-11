module CraftsJuMPExt

using Crafts, JuMP, SparseArrays

# Snap a float direction to a small exact rational; the true generators of count cones have small
# denominators, so this usually recovers them exactly (but nothing certifies that).
_snap(v) = rationalize.(BigInt, v ./ maximum(abs, v); tol=1e-8)

# Float support-LP oracle over the cone `{x : A*x <= 0}` (linearity rows tight). The model is
# built once; each probe only swaps the objective, so warm starts carry over between solves.
#
# The probed directions must be exact facet candidates (they are: projectcone probes facets of
# the exact hull of the current generators): over a cone the support value is 0 or +∞,
# discontinuous in the direction, so a direction tilted off a facet by float noise would read as
# unbounded no matter the tolerance.
function Crafts._lporacle(optimizer::Union{Type,Function,MOI.OptimizerWithAttributes},
                          P, A, linearity)
    m, n = size(A)
    islin = falses(m)
    islin[linearity] .= true
    Af = sparse(Float64.(A))
    Pf = Float64.(P)

    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[1:n])
    @constraint(model, Af[.!islin, :] * x .<= 0)
    @constraint(model, Af[islin, :] * x .== 0)

    return function (c)
        w = Float64.(transpose(P) * c)
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

end
