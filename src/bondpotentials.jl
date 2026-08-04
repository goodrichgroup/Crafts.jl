"""
    harmonicpotential(As::AbstractMatrix, Bs::AbstractMatrix; k)

Construct a harmonic potential with patch positions determined by matrices `A` and `B`,
and spring constant `k`.
"""
function harmonicpotential(A::AbstractMatrix, B::AbstractMatrix; k)
    n = size(A, 2)
    K = k / n

    function energy(p1::Pose, p2::Pose)
        aᵢ = (p1 * a for a in eachcol(A))
        bᵢ = (p2 * b for b in eachcol(B))
        return 0.5 * K * sum(sum((a[i] - b[i])^2 for i in eachindex(a,b)) for (a, b) in zip(aᵢ, bᵢ))
    end

    return energy
end

"""
    harmonicpotential(As::AbstractMatrix, Δp::Pose; k::Real)

Construct a harmonic potential with patch positions determined by matrices `A` and a equilibrium pose `Δp`, with spring
constants `k`.
"""
function harmonicpotential(As::AbstractMatrix, Δp::Pose; k::Real)
    Bs = reduce(hcat, inv(Δp) * a for a in eachcol(As))
    return make_harmonicpotential(As, Bs; k)
end

"""
    harmonicpotential(σ::Real; k::Real)

Construct a harmonic potential with uniform patch positions along a line of length `σ` and total spring constant `k`.
"""
function harmonicpotential(σ::Real=1; k::Real)
    t = σ^2 / 12
    A = [0 0; 0 t]

    function energy(p1::Pose, p2::Pose)
        p2 = p2 * Roly.standard_offset(Pose{2}()) # TODO: make nicer
        Δx = p1.psi' * (p2.x - p1.x)
        ΔR = p1.psi' * p2.psi
       
        return k * (t + 0.5 * sum(Δx[i]^2 for i in eachindex(Δx)) - tr(ΔR * A))
    end

    return energy
end

"""
    harmonicpotential(σx::Real, σy::Real; k::Real)

Construct a harmonic potential with uniform patch positions along a sheet of dimensions `σx` and `σy` and total spring
constant `k`.
"""
function harmonicpotential(σx::Real, σy::Real; k::Real)
    A = 1 / 12 * [0 0 0; 0 σx^2 0; 0 0 σy^2]
    t = tr(A)

    function energy(p1::Pose, p2::Pose)
        p2 = p2 * Roly.standard_offset(Pose{3}()) # TODO: make nicer
        Δx = p1.psi' * (p2.x - p1.x)
        ΔR = p1.psi' * p2.psi
       
        return k * (t + 0.5 * sum(Δx[i]^2 for i in eachindex(Δx)) - tr(ΔR * A))
    end
    return energy
end

abstract type BondPotential end

struct RigidSpringPotential{D,F} <: BondPotential
    Sij::Matrix{SMatrix{D,D,F}}
    abar::Vector{SVector{D,F}}
    k::Matrix{F}
    color_independent::Bool
end

function (rsp::RigidSpringPotential{D})(p1::Pose{D}, p2::Pose{D}; color1, color2) where {D}
    if rsp.color_independent
        color1 = color2 = 1
    end

    S11 = rsp.Sij[color1, color1]
    S22 = rsp.Sij[color2, color2]
    S12 = rsp.Sij[color1, color2]
    a1 = p1.psi * abar[color1]
    a2 = p2.psi * abar[color2]
    k = rsp.k[color1, color2]

    Δx = p2.x - p1.x
    ΔR = p1.psi' * p2.psi

    return k/2 * (tr(S11) + tr(S22) + normsq(Δx - (a1 - a2)) - 2*tr(S12*ΔR))
end