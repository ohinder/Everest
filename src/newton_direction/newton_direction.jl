abstract abstract_newton_direction;
abstract abstract_residuals;

type class_theta
    mu::Float64
    dual::Float64
    primal::Float64
end


include("class_homogeneous_residuals.jl")
include("class_homogeneous_newton.jl")


