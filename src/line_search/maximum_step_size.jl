
#######################
## MAXIMUM STEP
#######################

# while maintaining 0.1 * μ <= s' x 


function maximum_step(vars::class_variables, direction::class_variables)
	try
		alpha = Inf;
		alpha = maximum_step(alpha, x(vars), x(direction));
		alpha = maximum_step(alpha, s(vars), s(direction));
		alpha = maximum_step(alpha, tau(vars), tau(direction));
		alpha = maximum_step(alpha, kappa(vars), kappa(direction));

		return alpha;
	catch e
		println("ERROR in maximum_step")
		throw(e)
	end
end

function maximum_step(alpha::Float64, array_point::Array{Float64}, array_direction::Array{Float64})
	try
		@assert(length(array_point) == length(array_direction))
		for i in 1:length(array_point)
			alpha = maximum_step(alpha, array_point[i], array_direction[i]);
		end

		return alpha
	catch e
		println("ERROR in maximum_step(::Float64,::Array,::Array)")
		throw(e)
	end
end

function maximum_step(alpha::Float64, var::Float64, dir::Float64)
	try
		candidate_alpha = -var/dir;
		if candidate_alpha >= 0
			alpha = min(alpha, candidate_alpha)
		end

		return alpha
	catch e
		println("ERROR in maximum_step(::Float64,::Float64,::Float64)")
		throw(e)
	end
end
