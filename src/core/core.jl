abstract internal_AbstractNLPEvaluator;

include("class_variables.jl")
include("class_nlp_cache.jl")
include("class_quadratic_program.jl")
#include("class_qcqp.jl")
include("class_nlp.jl")
include("class_nlo.jl")
include("class_settings.jl")


function validate_dimensions( nlp_eval::internal_AbstractNLPEvaluator, vars::class_variables )
	try
    @assert(n(vars) == n(nlp_eval))
    @assert(m(vars) == m(nlp_eval))
 	catch e
	  println("ERROR validate_dimensions")
		throw(e)
	end
end

function validate_dimensions( nlp_val::class_nlp_cache )
    try
        num_vars = nlp_val.n
        num_constraints = nlp_val.m

        @assert( size( nlp_val.val_hesslag_prod ) == (num_vars,  num_vars) )
        @assert( size( nlp_val.val_jac_a ) == ( num_constraints,  num_vars) )
        @assert( size( nlp_val.val_a ) == ( num_constraints ,))
        @assert( size( nlp_val.val_gradlag ) == ( num_vars,))
        @assert ( size( nlp_val.val_gradc ) == ( num_vars,) )
        @assert ( size( nlp_val.val_b ) == ( num_constraints ,) )
    catch e
        println("ERROR validate_dimensions")
        throw(e)
	  end
end
