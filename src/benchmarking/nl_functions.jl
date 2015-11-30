function sum_sqrt(x::Array{Float64,1})
    return sum(sqrt(x+1))
end

function d_sum_sqrt(x::Array{Float64,1})
    return 0.5 * (x+1).^(-0.5)
end

function dd_sum_sqrt(x::Array{Float64,1})
    return spdiagm(-0.25 * (x+1).^(-1.5))
end

function nl_sqrt()
    return class_nl_function(sum_sqrt,d_sum_sqrt,dd_sum_sqrt);
end
