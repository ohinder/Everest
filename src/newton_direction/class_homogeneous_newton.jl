type class_homogeneous_newton <: abstract_newton_direction
	delta::Float64
	direction::class_variables
	residuals::class_homogeneous_residuals
  nlp_vals::class_nlp_cache

	K::SparseMatrixCSC{Float64,Int64}
  W::woodbury_identity
  W_updated::Bool
  F_updated::Bool

	linear_system_solver::abstract_linear_system_solver

  function class_homogeneous_newton()
      return new();
  end
end

###################
# methods
###################

function initialize_newton!(newt::class_homogeneous_newton, nlp_eval::internal_AbstractNLPEvaluator, vars::class_variables, settings::class_settings)
      try
				dim = n(vars) + m(vars) + 1;
				newt.K = spzeros(dim, dim);

				newt.direction = class_variables(n(vars),m(vars));

				newt.linear_system_solver = settings.linear_system_solver;
				initialize!(newt.linear_system_solver, newt.K);

        newt.nlp_vals = class_nlp_cache();

				newt.residuals = class_homogeneous_residuals();
        update_residuals!(newt.residuals, nlp_eval, vars, newt);
        newt.W_updated = false;

        newt.delta = 0.0;
			catch e
				println("ERROR in class_homogeneous_newton.intialize")
				throw(e)
			end
end


function update_newton!(newt::class_homogeneous_newton, vars::class_variables, settings::class_settings)
    try
      nlp_vals = newt.nlp_vals;

      start_advanced_timer("Matrix creation");

      val_x_scaled = x_scaled(vars);
      H = nlp_vals.val_hesslag_prod;
      h = sparse(H * val_x_scaled);
      c = nlp_vals.val_gradc;
      A = nlp_vals.val_jac_a;
      b = sparse(nlp_vals.val_b);

      #res.H += newt.delta * speye(length(val_x_scaled))
      D_x = H + spdiagm( s(vars) ./ x(vars) ) #+ ( newt.delta + settings.diagonal_modification ) * speye( n(vars) )
      D_g = sparse(val_x_scaled' * H * val_x_scaled + kappa(vars) / tau(vars)) #+ settings.diagonal_modification;
      #D_g = val_x_scaled' * nlp_vals.val_hesslag_prod * val_x_scaled + vars.kappa() / vars.tau() + this.delta;
      D_z = settings.diagonal_modification * speye(m(vars));
      #v_1 = sparse(-c - h);
      #v_2 = sparse(c - h);
      #v_3 = a - A * val_x_scaled;

     # newt.K_true = [
     #   [ D_x  	c - h	-A' 	];
     #   [ -c -h' 	D_g 	b' 	];
     #   [ A 	 -b	    D_z	]
     #   ];

      newt.K[:,:] = [
        [ D_x  		-h 	   A' 	];
        [ -h'	    D_g 	 -b' 	];
        [ A 		-b	   -D_z	]
        ];

      newt.F_updated = false

      pause_advanced_timer("Matrix creation");
  catch e
      println("ERROR in class_homogeneous_newton.update_newton!")
      throw(e)
  end
end

function update_newton_diag!(newt::class_homogeneous_newton, vars::class_variables, settings::class_settings)
    H = newt.nlp_vals.val_hesslag_prod;
    D_x_diag = diag(H) + s(vars) ./ x(vars) + newt.delta

    val_x_scaled = x_scaled(vars);
    D_g = val_x_scaled' * H * val_x_scaled + kappa(vars) / tau(vars) + newt.delta * norm(val_x_scaled,2)^2;
    diag_mod = [ D_x_diag; D_g; -settings.diagonal_modification*ones(m(vars))];

    for i = 1:size(newt.K,1)
        newt.K[i,i] = diag_mod[i];
    end

    return factorize_newton!(newt, vars)
end

function update_newton_diag_affine!(newt::class_homogeneous_newton, vars::class_variables, settings::class_settings)
    H = newt.nlp_vals.val_hesslag_prod;
    D_x_diag = diag(H) + s(vars) ./ x(vars) + newt.delta * x(vars).^(-2)

    val_x_scaled = x_scaled(vars);
    D_g = val_x_scaled' * H * val_x_scaled + kappa(vars) / tau(vars) + newt.delta / tau(vars)^2 #+ newt.delta * norm(val_x_scaled,2)^2;
    diag_mod = [ D_x_diag; D_g; -settings.diagonal_modification*ones(m(vars))];

    for i = 1:size(newt.K,1)
        newt.K[i,i] = diag_mod[i];
    end

    return factorize_newton!(newt, vars)
end

function factorize_newton!(newt::class_homogeneous_newton, vars::class_variables)
      start_advanced_timer("Factor");
      inertia = ls_factor(newt.linear_system_solver, n(vars) + 1, m(vars));
      pause_advanced_timer("Factor");

      newt.W_updated = false;
      newt.F_updated = true;

      return inertia
end

function form_woodbury!(newt::class_homogeneous_newton, vars::class_variables)
    try
        @assert(newt.F_updated)

        start_advanced_timer("woodbury_factor");
        c = newt.nlp_vals.val_gradc;

        vec = [c; zeros( 1 + m(vars) )];
        len = n(vars) + m(vars) + 1;
        U = [-e_( 1 + n(vars), len) vec];
        V = [vec e_( 1 + n(vars), len)]';

        newt.W = woodbury_identity(newt.linear_system_solver, U, V);

        newt.W_updated = true;
        pause_advanced_timer("woodbury_factor");
  catch e
      println("ERROR form_woodbury")
      throw(e)
  end
end

function compute_newton_direction!(newt::class_homogeneous_newton, vars::class_variables, theta::class_theta)
			try
        num_vars = n(vars);
        num_const = m(vars);
        res = newt.residuals;

				mu = res.mu;

				xs = -x(vars) .* s(vars) + theta.mu * mu * ones(num_vars);
				tk = -tau(vars) * kappa(vars) + theta.mu * mu;

				rhs =
        [	(1.0 - theta.dual) * res.r_D + xs ./ x(vars);
				(1.0 - theta.dual) * res.r_G + tk / tau(vars);
				(1.0 - theta.primal) * res.r_P ];

				#linear_system_solver.ls_solve!(rhs,sol);
        #start_advanced_timer("Factor2");
				#fac = lufact(newt.K_true)
				#pause_advanced_timer("Factor2");

        @assert(newt.W_updated == true)
        start_advanced_timer("Solve");
        #sol = (newt.K + U * newt.W.V) \ rhs;

        tol = 1e-6
        sol = false;
        try
          sol = ls_solve(newt.W, rhs);
          err = norm(evaluate(newt.W, sol) - rhs,1)/norm(rhs,1);
          @assert(err < tol)
        catch e
            warn("numerical instability using Woodbury, using direct factorization instead of woodbury.")
            sol = ls_solve_direct(newt.W, rhs)
            err = norm(evaluate(newt.W, sol) - rhs, 1)/norm(rhs,1);
            if (err > tol)
                warn("numerical stability using direct factorization, computation skipped")
                return false
            end
        end
        pause_advanced_timer("Solve");

        #println(full(newt.K_true))
        #abs_error = norm(newt.K_true * sol - rhs,1)
        #println("absolute_error=",abs_error)

        # update direction from linear system solution
				dir = newt.direction;
				x(dir, sol[1:num_vars] );
				tau(dir, sol[ num_vars + 1] );
				y(dir, -sol[( num_vars + 2):( num_const +  num_vars + 1)] );

				s(dir,(xs - s(vars) .* x(dir)) ./ x(vars) );
				kappa(dir, (tk - kappa(vars) * tau(dir)) / tau(vars) );

				# is this direction valid ?
				check_for_wrong_vals(dir);

        return true
			catch e
				println("ERROR in class_newton_solver.compute_newton_direction!")
        #@show full(newt.K + newt.W.U * newt.W.V)
				throw(e)
			end
		end

function mu(newton_direction::class_homogeneous_newton, vars::class_variables)
    return (dot(s(vars),x(vars)) + dot(tau(vars),kappa(vars)))/(n(vars) + 1)
end

