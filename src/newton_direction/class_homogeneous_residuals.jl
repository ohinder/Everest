type class_homogeneous_residuals <: abstract_residuals
	r_G::Float64
	r_D::Array{Float64,1}
	r_P::Array{Float64,1}

	r_norm::Float64
	r_D_norm::Float64
	r_G_norm::Float64
	r_P_norm::Float64

	r_D_norm_scaled::Float64
	r_G_norm_scaled::Float64
	r_P_norm_scaled::Float64

	mu::Float64
	scaled_mu::Float64

  val_c::Float64

	primal_norm::Float64
	dual_norm::Float64

	primal_infeas_sign::Int64
	dual_infeas_sign::Int64
	primal_infeas_norm::Float64
	dual_infeas_norm::Float64

  function class_homogeneous_residuals()
      return new();
  end
end

function finite_difference_merit_function(res::class_homogeneous_residuals, theta::class_theta, nlp, newton_solver, intial_point::class_variables, direction::class_variables)
		new_vars = deepcopy(intial_point)

		alpha = 1e-4;
		move!(new_vars, alpha, direction);
		orginal_mf = current_merit_function(res)
		val_mf = proximal_merit_function!(deepcopy(res), nlp, new_vars, newton_solver, intial_point);

		finite_dif = (val_mf - orginal_mf)/alpha

		computed_dif = merit_function_derivative(res, theta)

		@show finite_dif, computed_dif
end


function update_residuals!(res::class_homogeneous_residuals, nlp_eval::internal_AbstractNLPEvaluator, vars::class_variables, newt::abstract_newton_direction)
		update_residuals!(res, nlp_eval, vars, newt, vars)
end

function update_residuals!(res::class_homogeneous_residuals, nlp_eval::internal_AbstractNLPEvaluator, vars::class_variables, newt::abstract_newton_direction, intial_point::class_variables)
			try
				start_advanced_timer("residuals");

        nlp_vals = newt.nlp_vals;
        update_nlp_cache!(nlp_vals, nlp_eval, vars)

				res.mu = mu(newt,vars);
        res.val_c = nlp_vals.val_c

				val_x_scaled = x_scaled(vars);
        val_s = s(vars);
        val_kappa = kappa(vars);
        val_tau = tau(vars);


				prox_mod_x = -newt.delta_mod[1:n(vars)] .* (x(vars) - x(intial_point))
				res.r_D = s(vars) - tau(vars) * nlp_vals.val_gradlag + prox_mod_x;

				prox_mod_tau = -newt.delta_mod[n(vars)+1] * (tau(vars) - tau(intial_point))
				res.r_G = kappa(vars) + dot(nlp_vals.val_gradlag, x(vars)) + dot(nlp_vals.val_a, y(vars)) + prox_mod_tau;

				prox_mod_y = - newt.delta_mod[(n(vars)+2):(n(vars) + 1 + m(vars))] .* (y(vars) - y(intial_point))
				res.r_P = -tau(vars) * nlp_vals.val_a + prox_mod_y;

				res.r_D_norm = norm(res.r_D,1);
				res.r_G_norm = abs(res.r_G);
				res.r_P_norm = norm(res.r_P,1);
				res.r_norm = res.r_D_norm + res.r_G_norm + res.r_P_norm

				# other stuff
				scale = kappa(vars) + tau(vars);
				res.r_D_norm_scaled = res.r_D_norm/scale;
				res.r_G_norm_scaled = res.r_G_norm/scale;
				res.r_P_norm_scaled = res.r_P_norm/scale;

				res.scaled_mu = res.mu/(kappa(vars) + tau(vars));

				res.primal_norm = norm(res.r_P,1)/tau(vars);
				res.dual_norm = norm(res.r_D,1)/tau(vars);

				primal_infeas_obj = dot(nlp_vals.val_b, y(vars));
				dual_infeas_obj = -dot(x_scaled(vars), nlp_vals.val_gradc);
				res.primal_infeas_sign = sign(primal_infeas_obj);
				res.dual_infeas_sign = sign(dual_infeas_obj);

				res.primal_infeas_norm = norm(s(vars) + nlp_vals.val_jac_a' * y(vars),1)/abs(primal_infeas_obj);
 				res.dual_infeas_norm = norm(nlp_vals.val_a,1)/abs(dual_infeas_obj);

				pause_advanced_timer("residuals");
        #println(res.dual_infeas_norm, norm(nlp_vals.val_a,1)/nlp_vals.val_c)
			catch e
				  println("ERROR in class_residuals.update")
				  throw(e)
			end
end
