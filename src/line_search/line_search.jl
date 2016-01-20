function line_search(newton_solver::abstract_newton_direction, nlp::internal_AbstractNLPEvaluator, intial_point::class_variables, direction::class_variables, theta::class_theta)
	try
		alpha_max = maximum_step(intial_point, direction);
		damper = 0.01 #mu(newton_solver, intial_point))
		alpha = min(1.0,(1 -  damper) * alpha_max);

		res = newton_solver.residuals

		merit_function_diff = merit_function_derivative(res, theta)
		#finite_difference_merit_function(res, theta, nlp, newton_solver, intial_point, direction)

		orginal_merit_func = current_merit_function(res)

		return_point = None

		for i = 1:100
				new_vars = deepcopy(intial_point)
				move!(new_vars, alpha, direction);

				# why do I need to deep copy the residuals ???????
				new_merit_func = proximal_merit_function!(deepcopy(res), nlp, new_vars, newton_solver, intial_point)
				actual_gain = orginal_merit_func - new_merit_func

				v = new_merit_func/orginal_merit_func

				if actual_gain > 1e-4 * merit_function_diff * alpha
						return new_vars, alpha, true, v
				elseif alpha < 1e-5
						return new_vars, alpha, false, v
				end

				alpha = 0.99 * alpha * max(alpha,0.01)
		end

		error("failure")
	catch e
		println("ERROR in line_search")
		throw(e)
	end
end

#function merit_function_derivative()
#
#end

#function merit_function(newt::class_homogeneous_newton)
#
#end

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
