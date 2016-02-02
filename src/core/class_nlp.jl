type class_nl_function
    n::Int64
    value::Function
    gradient::Function
    hessian::Function

    function class_nl_function(num_var::Int64,val::Function,grad::Function,hess::Function)
        return new(num_var,val,grad,hess)
    end
end

type class_nlp <: internal_AbstractNLPEvaluator
  # linear program with non-linear objective
  # min f(x)
  # g_i(x) = 0
  # x >= 0

	_n::Int64 # number of variables
	_m::Int64 # number of constraints
  obj::class_nl_function
  constraints::Array{class_nl_function,1}

	function class_nlp(number_variables::Int64, obj::class_nl_function, constraints::Array{class_nl_function,1})
		return new(number_variables, length(constraints), obj, constraints);
	end
end

################
# METHODS
################

# evaluate objective
function internal_eval_c(nlp::class_nlp, x::Array{Float64,1})
    try
      return nlp.obj.value(x);
    catch e
        println("ERROR internal_eval_c")
        throw(e)
    end
end

# evalutate constraints
function internal_eval_a(nlp::class_nlp, x::Array{Float64,1})
    try
        num_constr = length(nlp.constraints)
        a = zeros(num_constr)
        for i = 1:num_constr
            a[i] = nlp.constraints[i].value(x)
        end

        return a;
    catch e
        println("ERROR internal_eval_a")
        throw(e)
    end
end

# evaluate gradient of constraints
function internal_eval_jac_a(nlp::class_nlp, x::Array{Float64,1}) # J
    try
      J = spzeros(m(nlp), n(nlp))

      for i = 1:m(nlp)
          J[i,:] = nlp.constraints[i].gradient(x)'
      end

      return J
    catch e
        println("ERROR internal_eval_jac_a")
        throw(e)
    end
end

# hessian of lagrangian
function internal_eval_hesslag_prod(nlp::class_nlp, x::Array{Float64,1}, y::Array{Float64,1})
    try
      H = nlp.obj.hessian(x)

      for i = 1:m(nlp)
          H = H - y[i] * nlp.constraints[i].hessian(x)
      end

      return H
    catch e
        println("ERROR internal_eval_hesslag_prod")
        throw(e)
    end
end

# gradient of lagrangian
function internal_eval_gradlag(nlp::class_nlp, x::Array{Float64,1}, y::Array{Float64,1})
    try
        g = internal_eval_gradc(nlp, x)

        for i = 1:m(nlp)
            g = g - y[i] * nlp.constraints[i].gradient(x)
        end
        return g
    catch e
        println("ERROR internal_eval_gradlag")
        throw(e)
    end
end

# gradient of f
function internal_eval_gradc(nlp::class_nlp, x::Array{Float64,1})
    try
        return nlp.obj.gradient(x);
    catch e
        println("ERROR internal_eval_gradc")
        throw(e)
    end
end

# (\nabla g * x - g)
function internal_eval_b(nlp::class_nlp, x::Array{Float64,1})
    try
        return internal_eval_jac_a(nlp, x) * x - internal_eval_a(nlp, x);
    catch e
        println("ERROR internal_eval_b")
        throw(e)
    end
end

function n(nlp::class_nlp)
    return nlp._n;
end

function m(nlp::class_nlp)
    return nlp._m;
end
