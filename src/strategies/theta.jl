function predictor_corrector(newton_solver::abstract_newton_direction, vars::class_variables, settings::class_settings)
	try
		# predictor
		gamma = 0.0;
		compute_newton_direction!(newton_solver, vars, class_theta(gamma,gamma,gamma)); # eta_P, eta_D, eta_G ???

		predictor_vars, = line_search(vars, newton_solver.direction);

		intial_mu = mu(newton_solver,vars);
		predictor_mu = mu(newton_solver,predictor_vars);

		# corrector
		reduction = predictor_mu/intial_mu;
		gamma = reduction^3 #(reduction)^3 #*0.8;
		compute_newton_direction!(newton_solver, vars, class_theta(gamma,gamma,gamma));
		vars, alpha = line_search(vars, newton_solver.direction);

		return vars, alpha, gamma
	catch e
		println("ERROR in predictor_corrector")
		throw(e)
	end
end

function simple_gamma_strategy(newton_solver::abstract_newton_direction, vars::class_variables, settings::class_settings)
	try
    η = 1.0;
    compute_newton_direction!(newton_solver, vars, class_theta(1.0,1 - η, 1 - η));
    vars, alpha = line_search(vars, newton_solver.direction);


		return vars, alpha, 1 - η
	catch e
		println("ERROR in simple_gamma_strategy")
		throw(e)
	end
end

function simple_gamma_strategy2(newton_solver::abstract_newton_direction, vars::class_variables, settings::class_settings)
	try
    alpha = 0.0
    η = 1.0
    i = 0
    MAX_IT = 10;

    for i = 1:MAX_IT
      compute_newton_direction!(newton_solver, vars, class_theta(1.0, 1 - η, 1 - η));
      vars, alpha = line_search(vars, newton_solver.direction);

      if alpha < 0.1
          η = η/1.5
      else
          break
      end
    end

    if i == MAX_IT
        println("too small eta choosen")
    end


		return vars, alpha, 1 - η
	catch e
		println("ERROR in simple_gamma_strategy")
		throw(e)
	end
end

function balancing_gamma_strategy(newton_solver::abstract_newton_direction, vars::class_variables, settings::class_settings)
	try
    alpha = 0.0
    r_P = newton_solver.residuals.r_P_norm;
    r_D = newton_solver.residuals.r_D_norm + newton_solver.residuals.r_G_norm;
    mu = newton_solver.residuals.mu;
    m = max(r_P,r_D,mu)
    η_P = r_P/m;
    η_D = r_D/m;
    η_mu = mu/m

    println(1 - η_mu,  " ", 1 - η_D, " ", 1 - η_P)
    compute_newton_direction!(newton_solver, vars, class_theta(1 - η_mu,1 - η_D, 1 - η_P));
    vars, alpha = line_search(vars, newton_solver.direction);

		return vars, alpha, 1.0
	catch e
		println("ERROR in simple_gamma_strategy")
		throw(e)
	end
end


function hybrid_mu_strategy(newton_solver::abstract_newton_direction, vars::class_variables, settings::class_settings, used_delta::Float64)
	try
		if used_delta > settings.delta_min
			vars, alpha, gamma = simple_gamma_strategy2(newton_solver, vars, settings)
      #vars, alpha, gamma = balancing_gamma_strategy(newton_solver, vars, settings)
		else
			vars, alpha, gamma = predictor_corrector(newton_solver, vars, settings)
		end
		return vars, alpha, gamma
	catch e
		println("ERROR in hybrid_mu_strategy")
		throw(e)
	end
end



function residuals_gamma_strategy(newton_solver::abstract_newton_direction, vars::class_variables, qp::class_quadratic_program)
	# try to balance residuals
	### not implemented
end
