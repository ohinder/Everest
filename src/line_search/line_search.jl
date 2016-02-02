include("maximum_step_size.jl")

function merit_function_derivative(res::class_homogeneous_residuals, theta::class_theta)
		return (theta.mu - 1.0) * res.mu  + (theta.dual - 1.0) * (res.r_G_norm + res.r_D_norm) + (theta.primal - 1.0) * res.r_P_norm;
end

function is_pos(arr::Array{Float64,1})
		for el in arr
				if el <= 0
						return false
				end
		end

		return true
end

function sx(vars::class_variables)
		sx_vec = [s(vars) .* x(vars); kappa(vars) * tau(vars)]
		return minimum(sx_vec)/mean(sx_vec) > 1e-2
end


function safe_guard(value::Float64, vars::class_variables)
		if is_pos(x(vars)) && is_pos(s(vars)) && tau(vars) > 0 && kappa(vars) > 0
			if sx(vars)
					return value
			else
					#print("SX!!")
					return Inf
			end
		else
			return Inf
		end
end

function current_merit_function(res::class_homogeneous_residuals, vars::class_variables)
		return safe_guard(res.mu + res.r_norm, vars)
end

function proximal_merit_function!(res::class_homogeneous_residuals, nlp_eval::internal_AbstractNLPEvaluator, vars::class_variables, newt::abstract_newton_direction, intial_point::class_variables)
		update_residuals!(res, nlp_eval, vars, newt, intial_point)
		return current_merit_function(res, vars)
end

function objective_merit_function(vars::class_variables, nlp::internal_AbstractNLPEvaluator, τ_k::Float64, μ_k::Float64)
		# (τ_k)^2 * c(𝑥/τ) - μ * ( sum(log(𝑥)) + log(τ) )
		obj = internal_eval_c(nlp, x_scaled(vars))
		value = (τ_k)^2 * obj - μ_k * ( sum(log(x(vars))) + log(tau(vars)) )
		return value
end

#function derivative_merit_function

function weighted_objective_constraints_merit_function(vars::class_variables, intial_point::class_variables, nlp::internal_AbstractNLPEvaluator, τ_k::Float64, μ_k::Float64)
		 y_mod = 1e-8 * (y(vars) - y(intial_point))
		 primal = norm(tau(vars) * internal_eval_a(nlp, x_scaled(vars)) + y_mod,1)
		 obj = objective_merit_function(vars, nlp, τ_k, μ_k)
		 #@show primal, obj
		 value = primal
		 return safe_guard(value,vars)
end

function step(newton_solver::abstract_newton_direction, nlp::internal_AbstractNLPEvaluator, intial_point::class_variables, direction::class_variables, theta::class_theta)
	try
		alpha_max = maximum_step(intial_point, direction);
		alpha = 0.5 * min(1.0, 0.8 * alpha_max);
		new_vars = deepcopy(intial_point)
		move!(new_vars, alpha, direction);

		res = newton_solver.residuals
		orginal_merit_func = proximal_merit_function!(res, nlp, intial_point, newton_solver, intial_point)
		new_merit_func = proximal_merit_function!(deepcopy(res), nlp, new_vars, newton_solver, intial_point)
		v = new_merit_func/orginal_merit_func

		return new_vars, alpha, true, v
	catch e
		println("ERROR in step")
		throw(e)
	end
end

function line_search(newton_solver::abstract_newton_direction, nlp::internal_AbstractNLPEvaluator, intial_point::class_variables, direction::class_variables, theta::class_theta)
	try
		alpha_max = maximum_step(intial_point, direction);
		alpha = min(1.0, 0.9 * alpha_max);

		res = newton_solver.residuals

		merit_function_diff = merit_function_derivative(res, theta)

		if merit_function_diff >= 0.0
				println(merit_function_diff)
				error("non-decreasing direction")
		end

		#finite_difference_merit_function(res, theta, nlp, newton_solver, intial_point, direction)

		#orginal_merit_func = current_merit_function(res)
		orginal_merit_func = proximal_merit_function!(res, nlp, intial_point, newton_solver, intial_point)

		actual_gain = 0.0
		new_merit_func = 0.0
		new_vars = intial_point
		for i = 1:20
				new_vars = deepcopy(intial_point)
				move!(new_vars, alpha, direction);

				# why do I need to deep copy the residuals ???????
				new_res = deepcopy(res)
				new_merit_func = proximal_merit_function!(new_res, nlp, new_vars, newton_solver, intial_point)

				# try non-linear update
				#
				nonlinear_vars = deepcopy(new_vars)
				non_linear_update(new_res, res, nonlinear_vars, (1 - alpha * (1 - theta.mu)))
				nonlinear_merit_func = proximal_merit_function!(new_res, nlp, nonlinear_vars, newton_solver, intial_point)

				is_non_linear_update = false
				if nonlinear_merit_func < new_merit_func
						new_vars = nonlinear_vars
						is_non_linear_update = true
						new_merit_func = nonlinear_merit_func
				end
				#

				# what is the gain?
				actual_gain = new_merit_func - orginal_merit_func
				v = new_merit_func/orginal_merit_func

				if actual_gain <= 1e-4 * merit_function_diff * alpha || alpha < 1e-4
					return new_vars, alpha, true, v
				end

				alpha = 0.5 * alpha
		end


		@show alpha
		@show actual_gain, merit_function_diff * alpha
		@show orginal_merit_func, new_merit_func
		@show norm(direction._v)
		@show norm(intial_point._v)
		@show norm(new_vars._v)
		@show norm(y(direction))
		@show norm(s(direction))
		@show tau(direction), kappa(direction)

		error("failure of line search")
	catch e
		println("ERROR in line_search")
		throw(e)
	end
end

function exact_line_search(newton_solver::abstract_newton_direction, nlp::internal_AbstractNLPEvaluator, intial_point::class_variables, direction::class_variables, theta::class_theta)
	try
		alpha_max = maximum_step(intial_point, direction);
		alpha = min(1.0, 0.9 * alpha_max);

		res = newton_solver.residuals

		merit_function_diff = merit_function_derivative(res, theta)

		if merit_function_diff >= 0.0
				error("non-decreasing direction")
		end

		#finite_difference_merit_function(res, theta, nlp, newton_solver, intial_point, direction)

		#orginal_merit_func = current_merit_function(res)
		orginal_merit_func = proximal_merit_function!(res, nlp, intial_point, newton_solver, intial_point)

		actual_gain = 0.0

		best_gain = actual_gain
		best_point = intial_point;
		best_alpha = 0.0
		best_v = NaN
		for i = 1:30
				new_vars = deepcopy(intial_point)
				move!(new_vars, alpha, direction);

				# why do I need to deep copy the residuals ???????
				new_res = deepcopy(res)
				new_merit_func = proximal_merit_function!(new_res, nlp, new_vars, newton_solver, intial_point)

				# try non-linear update
				#
					nonlinear_vars = deepcopy(new_vars)
					non_linear_update(new_res, res, nonlinear_vars, (1 - alpha * (1 - theta.mu)))
					nonlinear_merit_func = proximal_merit_function!(new_res, nlp, nonlinear_vars, newton_solver, intial_point)

					is_non_linear_update = false
					if nonlinear_merit_func < new_merit_func
							new_vars = nonlinear_vars
							is_non_linear_update = true
							new_merit_func = nonlinear_merit_func
					end
				#

				# what is the gain?
				actual_gain = new_merit_func - orginal_merit_func
				v = new_merit_func/orginal_merit_func

				if actual_gain <= best_gain
					best_gain = actual_gain
					best_point = new_vars;
					best_alpha = alpha
					best_v = v
				end

				alpha = 0.6 * alpha
		end
		if best_gain >= 0.0
				@show alpha
				@show actual_gain, merit_function_diff * alpha
				error("failure of line search")
		end

		return best_point, best_alpha, true, best_v
	catch e
		println("ERROR in line_search")
		throw(e)
	end
end


####################
#	TO DO!!!!
####################
function weighted_objective_constraint_line_search(newton_solver::abstract_newton_direction, nlp::internal_AbstractNLPEvaluator, intial_point::class_variables, direction::class_variables, theta::class_theta)
	try
		alpha_max = maximum_step(intial_point, direction);
		alpha = min(1.0, 0.9 * alpha_max);

		res = newton_solver.residuals

		merit_function_diff = -(1.0 - theta.primal) * res.r_P_norm

		if merit_function_diff >= 0.0
				@show merit_function_diff
				error("non-decreasing direction")
		end

		#finite_difference_merit_function(res, theta, nlp, newton_solver, intial_point, direction)

		τ_k = tau(intial_point)
		μ_k = theta.mu * mu(intial_point)
		orginal_merit_func = weighted_objective_constraints_merit_function(intial_point, intial_point, nlp, τ_k, μ_k)

		new_vars = deepcopy(intial_point)
		move!(new_vars, 1e-8, direction);

		inc_merit_func = weighted_objective_constraints_merit_function(new_vars, intial_point, nlp, τ_k, μ_k)

		@show merit_function_diff, (inc_merit_func - orginal_merit_func)/1e-8
		@show -μ_k * (sum(x(direction) ./ x(intial_point)) + tau(direction)/tau(intial_point))# + sum(1 ./ s(direction)) + 1./kappa(direction))
		@show dot(y(direction), res.r_P)
		#@show (1 - theta.mu) * dot(-y(direction), res.r_P) + dot(x(direction), res.r_D) + tau(direction) * res.r_G
		actual_gain = 1.0

		best_gain = actual_gain
		best_point = intial_point;
		best_alpha = 0.0
		best_v = NaN
		for i = 1:4
				new_vars = deepcopy(intial_point)
				move!(new_vars, alpha, direction);

				# why do I need to deep copy the residuals ???????
				new_res = deepcopy(res)
				new_merit_func = weighted_objective_constraints_merit_function(new_vars, intial_point, nlp, τ_k, μ_k)

				# what is the gain?
				actual_gain = new_merit_func - orginal_merit_func
				v = new_merit_func/orginal_merit_func

				if actual_gain <= best_gain
					best_gain = actual_gain
					best_point = new_vars;
					best_alpha = alpha
					best_v = v
				end

				alpha = 0.6 * alpha
		end
		if best_gain >= 1.0
				@show alpha
				@show actual_gain #, merit_function_diff * alpha
				@show norm(x(direction))
				@show norm(y(direction))
				@show norm(s(direction))
				@show tau(direction), kappa(direction)
				error("failure of weighted line search")
		end

		return best_point, best_alpha, true, best_v
	catch e
		println("ERROR in weighted_objective_constraint_line_search")
		throw(e)
	end
end



function oli_line_search(newton_solver::abstract_newton_direction, nlp::internal_AbstractNLPEvaluator, intial_point::class_variables, direction::class_variables, theta::class_theta)
	try
		# make sufficient progress on constraints
	catch

	end
end

#non_linear_update(previous_res, newton_solver.nlp_vals, deepcopy(new_vars), 0.99)
function non_linear_update(new_res::class_homogeneous_residuals, previous_res::class_homogeneous_residuals, point::class_variables, reduction_factor::Float64)
		s( point, reduction_factor * previous_res.r_D + (s(point) - new_res.r_D) )
		kappa( point, reduction_factor * previous_res.r_G + (kappa(point) - new_res.r_G) )
end

function s_heuristic()

end

#
