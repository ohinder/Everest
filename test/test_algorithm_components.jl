# TO DO:
# test termination criterion

using Base.Test
using FactCheck

function qp_test()
    Q = speye(3);
    A = sparse(ones(1,3));
    b = 3*ones(1);
    c = ones(3);
    return class_quadratic_program(A, b, c, Q), A, b, c, Q;
end

facts("core") do
  ###############################
  # test class_variables
  ###############################
  context("variables") do
    var = class_variables(3,3);

    x_values = [-1, 2, 2];
    var._v[var._x_ind] = x_values;

    @fact x(var) --> x_values
    @fact_throws ErrorException check_positive(var)

    x_values = [1, NaN, 2];
    var._v[var._x_ind] = x_values;
    @fact_throws ErrorException check_for_wrong_vals(var)

    x_values = [1, Inf, 2];
    var._v[var._x_ind] = x_values;
    @fact_throws ErrorException check_for_wrong_vals(var)

    x_values = [1, -Inf, 2];
    var._v[var._x_ind] = x_values;
    @fact_throws ErrorException check_for_wrong_vals(var)

    var = class_variables(3,3);
    @fact_throws ErrorException move!(var,-2.0,var)
  end

  qp, A, b, c, Q = qp_test();

  # test quadratic program and validate dimensions
  context("qp and validate dimensions") do
      x_var = [1.0,1.0,1.0]
      y_var = [1.0];

      @fact 4.5 --> internal_eval_c(qp, x_var)
      @fact [0.0] --> internal_eval_a(qp, x_var)
      @fact A --> internal_eval_jac_a(qp, x_var)
      @fact Q --> internal_eval_hesslag_prod(qp, x_var, y_var)
      @fact c + Q * x_var - A' * y_var --> internal_eval_gradlag(qp, x_var, y_var)

      vars2 = class_variables(3,3);
      @fact_throws ErrorException validate_dimensions(qp, vars2)
  end
end


begin
    qp, A, b, c, Q = qp_test();
    nlp_vals = class_nlp_cache();
    settings = class_settings();
    vars = class_variables(3,1);
    update_nlp_cache!(nlp_vals, qp, vars)
    validate_dimensions(qp, vars)

    newt = class_homogeneous_newton();
    # test newton direction
    facts("newton direction") do

      start_advanced_timer()
      begin
          newt = class_homogeneous_newton();
          initialize_newton!(newt, qp, vars, settings)
          res = class_homogeneous_residuals();
      end


      begin
          @fact 1.0 --> mu(newt, vars);

          validate_dimensions(qp, vars)
          initialize_newton!(newt, qp, vars, settings)
          update_newton!(newt, vars, settings)
          update_newton_diag!(newt, vars, settings)
          validate_dimensions( newt.nlp_vals )
          form_woodbury!(newt, vars)
          compute_newton_direction!(newt, vars, class_theta(0.5,0.5,0.5))
      end

  end

  facts("line search") do
    # test line search
    begin
        vars2= class_variables(5,3);
        direction = class_variables(5,3);
        x(direction, -x(vars2))
        vars2, alpha = line_search(vars2, direction);
        @fact alpha --> less_than(1)
    end
  end

  facts("strategies") do
      # test stratagies

      # test delta heuristics
      begin
          #println("TEST IPOPT INERTIA CORRECTION")
          validate_dimensions(qp, vars)
          delta, number_of_factors = ipopt_style_inertia_correction!(newt, vars, settings)
          @fact delta --> 0.0
          @fact number_of_factors --> 1
      end

       # test theta heursitics
      begin
          predictor_corrector(newt, vars, settings)
          simple_gamma_strategy(newt, vars, settings)
          hybrid_mu_strategy(newt, vars, settings, settings.delta_min + 1.0)
          hybrid_mu_strategy(newt, vars, settings, settings.delta_min)
      end

      pause_advanced_timer()

  end

    #begin
    #    vars = class_variables(3,1);
    #    A = sparse(ones(1,3));
    #    b = 3*ones(1);
    #    c = [2.0, 0.0, 1.0];
    #    qp = class_quadratic_program(A, b, c);
    #    status, vars = homogeneous_algorithm(qp, vars, settings)
    #    @fact norm(x_scaled(vars) - [0.0,3.0,0.0]) --> less_than(1e-6)
    #    @fact tau(vars) --> greater_than(0.1)
    #    @fact 0.0 --> less_than(kappa(vars))
    #    @fact kappa(vars) --> less_than(1e-6)
    #end


end
