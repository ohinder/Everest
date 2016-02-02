# step or line_search
#LINE_SEARCH_METHOD = step;
LINE_SEARCH_METHOD = line_search;
#LINE_SEARCH_METHOD = exact_line_search;
#LINE_SEARCH_METHOD = weighted_objective_constraint_line_search;


function predictor_corrector(newton_solver::abstract_newton_direction, nlp::internal_AbstractNLPEvaluator, vars::class_variables, settings::class_settings)
	try
		# predictor
		gamma = 0.0;
		theta = class_theta(gamma,gamma,gamma)
		compute_newton_direction!(newton_solver, vars, theta);

		predictor_vars, alpha, success, v = LINE_SEARCH_METHOD(newton_solver, nlp, vars, newton_solver.direction, theta);

		intial_mu = mu(newton_solver,vars);
		predictor_mu = mu(newton_solver,predictor_vars);

		# corrector
		reduction = predictor_mu/intial_mu;
		#gamma = 0.5
		gamma = maximum([0.01,reduction * minimum([reduction,0.9])])
		#gamma = v * minimum([v,0.75])

		res = newton_solver.residuals
		#if res.r_norm > res.mu * 10
		theta = class_theta(gamma,gamma,gamma)
		compute_newton_direction!(newton_solver, vars, theta);
		vars, alpha, = LINE_SEARCH_METHOD(newton_solver, nlp, vars, newton_solver.direction, theta);

		return vars, alpha, gamma
	catch e
		println("ERROR in predictor_corrector")
		throw(e)
	end
end

function simple_gamma_strategy(newton_solver::abstract_newton_direction, nlp::internal_AbstractNLPEvaluator, vars::class_variables, settings::class_settings)
	try
    η = 1.0;
		theta = class_theta(1.0,1 - η, 1 - η)
    compute_newton_direction!(newton_solver, vars, theta);
    vars, alpha = LINE_SEARCH_METHOD(newton_solver, nlp, vars, newton_solver.direction, theta);


		return vars, alpha, 1 - η
	catch e
		println("ERROR in simple_gamma_strategy")
		throw(e)
	end
end

function simple_gamma_strategy2(newton_solver::abstract_newton_direction, nlp::internal_AbstractNLPEvaluator, vars::class_variables, settings::class_settings)
	try
    alpha = 0.0
    η = 1.0
    i = 0
    MAX_IT = 3;

    for i = 1:MAX_IT
			theta = class_theta(1.0,1 - η, 1 - η)
      compute_newton_direction!(newton_solver, vars, theta);
      vars, alpha = LINE_SEARCH_METHOD(newton_solver, nlp, vars, newton_solver.direction, theta);
			#vars, alpha = step(newton_solver, nlp, vars, newton_solver.direction, theta);
			#vars, alpha = weighted_objective_constraint_line_search(newton_solver, nlp, vars, newton_solver.direction, theta);

      if alpha < 0.1
          η = 0.5 * η;
      else
          break
      end
    end

    if i == MAX_IT
        println("too small eta choosen")
    end


		return vars, alpha, 1 - η
	catch e
		println("ERROR in simple_gamma_strategy2")
		throw(e)
	end
end

function balancing_gamma_strategy(newton_solver::abstract_newton_direction, nlp::internal_AbstractNLPEvaluator, vars::class_variables, settings::class_settings)
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
		theta = class_theta(1 - η_mu,1 - η_D, 1 - η_P)
    compute_newton_direction!(newton_solver, vars, theta);
    vars, alpha = LINE_SEARCH_METHOD(newton_solver, nlp, vars, newton_solver.direction, theta);

		return vars, alpha, 1.0
	catch e
		println("ERROR in simple_gamma_strategy")
		throw(e)
	end
end


function hybrid_mu_strategy(newton_solver::abstract_newton_direction, nlp::internal_AbstractNLPEvaluator, vars::class_variables, settings::class_settings, used_delta::Float64)
	try
		if used_delta > settings.delta_min
			vars, alpha, gamma = simple_gamma_strategy2(newton_solver, nlp, vars, settings)
      #vars, alpha, gamma = balancing_gamma_strategy(newton_solver, vars, settings)
			return vars, alpha, gamma
		else
			vars, alpha, gamma = predictor_corrector(newton_solver, nlp, vars, settings)
			return vars, alpha, gamma
		end
	catch e
		println("ERROR in hybrid_mu_strategy")
		throw(e)
	end
end



function residuals_gamma_strategy(newton_solver::abstract_newton_direction, vars::class_variables, qp::class_quadratic_program)
	# try to balance residuals
	### not implemented
end
