type class_nl_function
    value::Function
    gradient::Function
    hessian::Function

    function class_nl_function(val::Function,grad::Function,hess::Function)
        return new(val,grad,hess)
    end
end

type class_nlo <: internal_AbstractNLPEvaluator
  # linear program with non-linear objective
  # min f(x)
  # A*x = b
  # x >= 0

	_n::Int64 # number of variables
	_m::Int64 # number of constraints
	_A::SparseMatrixCSC{Float64,Int64}
	_b::Array{Float64,1}
  obj::class_nl_function

	function class_nlo(A::SparseMatrixCSC{Float64,Int64},b::Array{Float64,1}, obj::class_nl_function)
		(m, n) = size(A)
		return new(n,m,A,b,obj);
	end
end

################
# METHODS
################

# evaluate objective
function internal_eval_c(nlo::class_nlo, x::Array{Float64,1})
    return nlo.obj.value(x);
end

# evalutate constraints
function internal_eval_a(nlo::class_nlo, x::Array{Float64,1})
    return nlo._A * x - nlo._b;
end

# evaluate gradient of constraints
function internal_eval_jac_a(nlo::class_nlo, x::Array{Float64,1}) # J
    return nlo._A;
end

# hessian of lagrangian
function internal_eval_hesslag_prod(nlo::class_nlo, x::Array{Float64,1}, y::Array{Float64,1})
    return nlo.obj.hessian(x)
end

# gradient of lagrangian
function internal_eval_gradlag(nlo::class_nlo, x::Array{Float64,1}, y::Array{Float64,1})
    return internal_eval_gradc(nlo, x) - nlo._A' * y;
end

# gradient of f
function internal_eval_gradc(nlo::class_nlo, x::Array{Float64,1})
    return nlo.obj.gradient(x);
end

# (\nabla g * x - g)
function internal_eval_b(nlo::class_nlo, x::Array{Float64,1})
    return nlo._b;
end

function n(nlo::class_nlo)
    return nlo._n;
end

function m(nlo::class_nlo)
    return nlo._m;
end
