abstract abstract_newton_direction;
abstract abstract_residuals;

type class_theta
    # mu = 1.0 <=> no mu change
    # mu = 0.0 <=> aggresively decrease mu
    mu::Float64
    dual::Float64
    primal::Float64
end


include("class_homogeneous_residuals.jl")
include("class_homogeneous_newton.jl")
