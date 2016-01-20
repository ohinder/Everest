# wrapper that rewrites the orginal problem as
# min f(x) + |x - x^k|_D^2
# to make it convex

# THIS NOT BEING USED CURRENTLY!!!

type class_proximal_eval <: internal_AbstractNLPEvaluator
    orginal_eval::internal_AbstractNLPEvaluator
    x::Array{Float64,1}
    D::SparseMatrixCSC{Float64,Int64} # the diagonal modification

    function class_proximal_eval(orginal_eval::internal_AbstractNLPEvaluator, starting_point::class_variables, D::SparseMatrixCSC{Float64,Int64})
        return new(orginal_eval, starting_point, D);
    end
end

# evaluate objective
function internal_eval_c(nlp::class_proximal_eval, x::Array{Float64,1})
    return internal_eval_c(nlp.orginal_eval,x) + 0.5 * (nlp.x - x) * D * (nlp.x - x);
end

# evalutate constraints
function internal_eval_a(nlp::class_proximal_eval, x::Array{Float64,1})
    return internal_eval_a(nlp.orginal_eval,x);
end

# evaluate gradient of constraints
function internal_eval_jac_a(nlp::class_proximal_eval, x::Array{Float64,1}) # J
    return internal_eval_jac_a(nlp.orginal_eval,x);
end

# hessian of lagrangian
function internal_eval_hesslag_prod(nlp::class_proximal_eval, x::Array{Float64,1}, y::Array{Float64,1})
    return internal_eval_hesslag_prod(nlp.orginal_eval,x,y) + D;
end

# gradient of lagrangian
function internal_eval_gradlag(nlp::class_proximal_eval, x::Array{Float64,1}, y::Array{Float64,1})
    return internal_eval_gradlag(nlp.orginal_eval,x,y) + D * (nlp.x - x);
end

# gradient of f
function internal_eval_gradc(nlp::class_proximal_eval, x::Array{Float64,1})
    return internal_eval_gradc(nlp.orginal_eval,x) + D * (nlp.x - x);
end

# (\nabla g * x - g)
function internal_eval_b(nlp::class_proximal_eval, x::Array{Float64,1})
    return internal_eval_b(nlp.orginal_eval,x)
end

function n(nlp::class_proximal_eval)
    return n(nlp.orginal_eval);
end

function m(nlp::class_proximal_eval)
    return m(nlp.orginal_eval)
end
