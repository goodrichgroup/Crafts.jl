using StaticArrays, Statistics, LinearAlgebra

function vertices2particleidxs(str, sys, vs)
    return unique!([Roly.vertex2particle(str, sys, v)[1] for v in vs])
end

function mobility(r)
    return inv(6π*r)
end
function rotneprager(Δx; ri=1, rj=1)
    R = norm(Δx)
    R2 = R^2
    P = Δx * Δx' / R2
    return (I + P + (ri^2 + rj^2) / R2 * (I/3 - P)) / (8π * R)
end
function frictiontensor(xs, rs)
    n = length(xs)
    V = 4/3 * π * sum(x->x^3, rs)
    B = zeros(3n, 3n)
    
    for i in 1:n
        B[1 + 3(i-1):3i, 1 + 3(i-1):3i] .= I(3) * mobility(rs[i])
        for j in i+1:n
            t = rotneprager(xs[j] - xs[i]; ri=rs[i], rj=rs[j])
            B[1 + 3(i-1):3i, 1 + 3(j-1):3j] .= t
            B[1 + 3(j-1):3j, 1 + 3(i-1):3i] .= t
        end
    end

    C = inv(B)

    U(i) = let (x,y,z) = (xs[i][1], xs[i][2], 0)
        [0 -z y;
         z 0 -x;
         -y x 0]
    end

    θt = sum(C[1 + 3(i-1):3i, 1 + 3(j-1):3j] for i in 1:n, j in 1:n)
    θtr = sum(U(i) * C[1 + 3(i-1):3i, 1 + 3(j-1):3j] for i in 1:n, j in 1:n)
    θrr_uncor = -sum(U(i) * C[1 + 3(i-1):3i, 1 + 3(j-1):3j] * U(j) for i in 1:n, j in 1:n)
    θrr = θrr_uncor + 6V * I

    θ = [θt θtr'; θtr θrr]
    return θ
end
function centerofdiffusion(xs, rs)
    D0 = inv(frictiontensor(xs, rs))

    dtr = @view D0[4:6, 1:3]
    drr = @view D0[4:6, 4:6]

    A = tr(drr) * I - drr
    b = [dtr[2,3] - dtr[3,2], dtr[3,1] - dtr[1,3], dtr[1,2] - dtr[2,1]]

    r_od = A \ b
    return r_od
end
function diffusiontensor(xs, rs; cod=nothing)
    if isnothing(cod)
        cod = centerofdiffusion(xs, rs)
    end
    D = inv(frictiontensor(xs .- Ref(cod), rs))
    return D
end
diffusionconstants(xs, rs; kwargs...) = let d = diffusiontensor(xs, rs; kwargs...)
    # Return Dlin, Drot
    (d[1, 1] + d[2, 2] + d[3, 3]) / 3, (d[4, 4] + d[5, 5] + d[6, 6]) / 3
end
function diffusionconstants(p::Polyform, sys::AssemblySystem; kwargs...)
    xs = [@SVector([x[1], x[2], 0]) for x in p.xs]
    rs = [sys.geometries[i].R_min for i in p.species]
    return diffusionconstants(xs, rs; kwargs...)
end

caparea(ϕ) = 4π * sin(ϕ / 2)^2
function brownianrate(; D, Dr1, δ1, Dr2, δ2, R)
    # for strongly non-convex bodies, it may be the case that R1 + R2 != R
    F1 = caparea(δ1) / 4π
    F2 = caparea(δ2) / 4π

    ξ1 = sqrt((1 + Dr1*R^2/D) / 2)
    ξ2 = sqrt((1 + Dr2*R^2/D) / 2)

    Λ1 = F1*(ξ1 + cot(δ1 / 2)) / (ξ1 + sin(δ1/2)*cos(δ1/2))
    Λ2 = F2*(ξ2 + cot(δ2 / 2)) / (ξ2 + sin(δ2/2)*cos(δ2/2))

    denom = Λ1*Λ2 + inv(inv(1 - Λ1)*inv(1 - Λ2) + inv(1 - Λ1) * inv(Λ2 - F2) + inv(1 - Λ2) * inv(Λ1 - F1))
    k = 4π*R*D*F1*F2 / denom 
    return k
end

function bindingrate(aggregate, sys, cut, components; Lsite=1)
    xs = [@SVector([x[1], x[2], 0]) for x in aggregate.xs]
    rs = [sys.geometries[i].R_min for i in aggregate.species]

    # Center of binding
    cob = [mean(Roly.get_sitepos(aggregate, sys, c.src) for c in cut); 0]

    # Site indices of the two components
    vs1, vs2 = components
    idxs1 = vertices2particleidxs(aggregate, sys, vs1)
    idxs2 = vertices2particleidxs(aggregate, sys, vs2)

    xs1, rs1 = xs[idxs1], rs[idxs1]
    xs2, rs2 = xs[idxs2], rs[idxs2]

    cod1 = centerofdiffusion(xs1, rs1)
    cod2 = centerofdiffusion(xs2, rs2)

    Dl1, Dr1 = diffusionconstants(xs1, rs1; cod=cod1)
    Dl2, Dr2 = diffusionconstants(xs2, rs2; cod=cod2)

    R = norm(cod2 - cod1)
    r1 = norm(cob - cod1)
    r2 = norm(cob - cod2)

    R1 = abs((r1^2 - r2^2 + R^2) / (2R))
    R2 = abs((r2^2 - r1^2 + R^2) / (2R))

    # Ω1 = Asite / R1^2
    # Ω2 = Asite / R2^2
    # δ1 = acos(1 - Ω1 / (2π))
    # δ2 = acos(1 - Ω2 / (2π))
    δ1 = atan(Lsite / (2R1))
    δ2 = atan(Lsite / (2R2))
    D = Dl1 + Dl2
    return brownianrate(; D, Dr1, δ1, Dr2, δ2, R)
end

function bindingrate_puretranslation(aggregate, sys, cut, components)
    xs = [@SVector([x[1], x[2], 0]) for x in aggregate.xs]
    rs = [sys.geometries[i].R_min for i in aggregate.species]

    # Site indices of the two components
    vs1, vs2 = components
    idxs1 = vertices2particleidxs(aggregate, sys, vs1)
    idxs2 = vertices2particleidxs(aggregate, sys, vs2)

    xs1, rs1 = xs[idxs1], rs[idxs1]
    xs2, rs2 = xs[idxs2], rs[idxs2]

    cod1 = centerofdiffusion(xs1, rs1)
    cod2 = centerofdiffusion(xs2, rs2)

    Dl1, _ = diffusionconstants(xs1, rs1; cod=cod1)
    Dl2, _ = diffusionconstants(xs2, rs2; cod=cod2)

    R = norm(cod2 - cod1)
    D = Dl1 + Dl2
    return 4π * D * R
end
