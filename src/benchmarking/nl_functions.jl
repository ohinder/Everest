# SQRT FUNCTION
# (x+1)^0.5
function sum_sqrt(x::Array{Float64,1})
    return sum(sqrt(x+1))
end

function d_sum_sqrt(x::Array{Float64,1})
    return 0.5 * (x+1).^(-0.5)
end

function dd_sum_sqrt(x::Array{Float64,1})
    return spdiagm(-0.25 * (x+1).^(-1.5))
end

function nl_sqrt(num_var::Int64)
    return class_nl_function(num_var,sum_sqrt,d_sum_sqrt,dd_sum_sqrt);
end

# LOG function
function sum_log(x::Array{Float64,1})
    return sum(log(x+1))
end

function d_sum_log(x::Array{Float64,1})
    return (x+1).^(-1.0)
end

function dd_sum_log(x::Array{Float64,1})
    return spdiagm(-(x+1).^(-2.0))
end

function nl_log(num_var::Int64)
    return class_nl_function(num_var::Int64,sum_log,d_sum_log,dd_sum_log);
end

# Quadratic
# c'*x + 0.5 * x*Q*x + b
function nl_quad(c::AbstractArray, Q::SparseMatrixCSC{Float64,Int64}, b::Float64)
    function value(x::Array{Float64,1})
        return (c' * x)[1] + 0.5 * (x' * Q * x)[1] + b
    end

    function gradient(x::Array{Float64,1})
        return (c + Q * x)[:]
    end

    function hessian(x::Array{Float64,1})
        return Q
    end

    return class_nl_function(length(c), value, gradient, hessian)
end

function sphere_function(center_point::Float64,radius::Float64)
    b = dot(center_point,center_point) - radius
    c = -2 * center_point

    return nl_quad(c, speye(length(c)), b)
end
