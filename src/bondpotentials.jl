function make_harmonicpotential(As::AbstractMatrix, Bs::AbstractMatrix; k::Real)
    n = size(As, 2)
    K = k / n

    function energy(p1::Pose, p2::Pose)
        aᵢ = (p1 * a for a in eachcol(As))
        bᵢ = (p2 * b for b in eachcol(Bs))
        return 0.5 * K * sum(sum((a[i] - b[i])^2 for i in eachindex(a,b)) for (a, b) in zip(aᵢ, bᵢ))
    end

    return energy
end

function make_harmonicpotential(As::AbstractMatrix, Δp::Pose; k::Real)
    Bs = reduce(hcat, inv(Δp) * a for a in eachcol(As))
    return make_harmonicpotential(As, Bs; k)
end

function make_harmonicpotential(σ::Real; k::Real)
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

function make_harmonicpotential(σ1::Real, σ2::Real; k::Real)
    A = 1 / 12 * [0 0 0; 0 σ1 0; 0 0 σ2]
    t = tr(A)

    function energy(p1::Pose, p2::Pose)
        p2 = p2 * Roly.standard_offset(Pose{2}()) # TODO: make nicer
        Δx = p1.psi' * (p2.x - p1.x)
        ΔR = p1.psi' * p2.psi
       
        return k * (t + 0.5 * sum(Δx[i]^2 for i in eachindex(Δx)) - tr(ΔR * A))
    end

    return energy
end