function ipopt_style_inertia_correction!(newt::abstract_newton_direction, vars::class_variables, settings::class_settings)
	try
    delta_style! = update_newton_diag!;

		MAX_IT = 25;
		j = 1;
		new_delta = 0.0;

		# try delta = 0.0
    old_delta = newt.delta;
    newt.delta = new_delta;
    update_newton!(newt, vars, settings);
    inertia = delta_style!(newt, vars, settings)

		if inertia == 1
			new_delta = 0.0; # no need to modify newton system
		else
			newt.delta = old_delta;

      if newt.delta <= settings.delta_min;
            newt.delta = settings.delta_start;
      end

			for j = 2:MAX_IT
				inertia = delta_style!(newt, vars, settings)

				if inertia == 1
            new_delta = newt.delta * settings.delta_decrease
            break
				elseif inertia == 0
            newt.delta = newt.delta * settings.delta_increase;
				elseif inertia == -1
					  error("numerical stability issues when computing KKT system !!!")
				else
					  error("inertia_corection")
				end
			end
		end

		if j == MAX_IT
			error("maximum iterations for inertia_corection reached")
		else
      form_woodbury!(newt, vars)
			return new_delta, j
		end


	catch e
		println("ERROR in ipopt_style_inertia_correction")
		throw(e)
	end
end


function iterative_trust_region!(newt::abstract_newton_direction, vars::class_variables, settings::class_settings)
    # Goal: choose δ such that Δ = |log(x) - log(x + d_x)|_∞ = 0.5
    # or |d_x/x|_∞ <= 0.5
    # or |d_x^2/x^2| <= 1.0
    #
    #
    # compute d'
    try
      newt.delta = 0.0;
      inertia = update_newton!(newt, vars, settings);

      if inertia
          best_δ = 0.0;
      else
          target = 2.0;
          best_val = Inf;
          δ_vals = 2.^linspace(0.0,4.0,20)

          for δ = δ_vals
              newt.delta = δ
              inertia = update_newton_diag!(newt, vars, settings);
              if inertia
                  form_woodbury!(newt, vars)
									
									try
		                  compute_newton_direction!(newt, vars, class_theta(0.99,0.5,0.5));
		                  alpha_max = maximum_step(vars, newt.direction);

		                  val = abs(log(alpha_max) - log(target));
		                  if best_val > val
		                      best_val = val;
		                      best_δ = δ;
		                  end
									catch e
											warn(e.msg)
											warn("computation skipped")
									end
              end
          end

          @assert(best_val < Inf)
      end

      form_woodbury!(newt, vars)
      return best_δ, 1
    catch e
        println("ERROR in delta.iterative_trust_region!")
        throw(e)
	  end
end

#vars = class_variables(1,1);
#dir = deepcopy(vars);
#x(dir,[-0.5])
#alpha_max = maximum_step(vars,dir)
