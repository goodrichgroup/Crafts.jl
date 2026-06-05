function n_species(M::AbstractMatrix)
    rows = sort(eachrow(M), by=sum)
    nμ = 0
    while sum(rows[nμ+1]) == 1
        nμ += 1
    end
    return nμ
end

normal_vec(x::SVector{2,F}) where F = SVector{2,F}(-x[2], x[1])

function infapprox(x, inf_val=99.9)
    return replace(x, Inf => inf_val, -Inf => -inf_val)
end