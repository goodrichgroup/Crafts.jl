

normal_vec(x::SVector{2,F}) where F = SVector{2,F}(-x[2], x[1])

function infapprox(x, inf_val=99.9)
    return replace(x, Inf => inf_val, -Inf => -inf_val)
end

normsq(v) = sum(x->x^2, v)