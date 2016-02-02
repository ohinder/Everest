
settings = class_settings();

context("nlp") do
    # - x - y

    a = sparse([-2.0, -2.0, 0.0])
    b = 0.0
    P = spdiagm([0.0, 0.0, 0.0])
    obj = nl_quad(a,P,b)

    # x^2 + y^2 - s - 4.0
    # x^2 + y^2 >= 4.0
    a1 = sparse([0.0, 0.0, -1.0])
    P1 = spdiagm([2.0, 2.0, 0.0])
    b1 = -4.0
    constr1 = nl_quad(a1,P1,b1)

    a1 = sparse([0.0, 0.0, -1.0])
    P1 = spdiagm([2.0, 2.5, 0.0])
    b1 = -4.0
    constr1_b = nl_quad(a1,P1,b1)

    # x + y + s - 1.0 = 0.0
    # x + y = 1.0
    a2 = sparse([1.0, 1.0, 0.0])
    P2 = spdiagm([0.0, 0.0, 0.0])
    b2 = -1.0 #-4.0
    constr2 = nl_quad(a2,P2,b2)

    # x + x^2 - y = 0
    # x + x^2 = y
    a3 = sparse([1.0, -1.0, 0.0])
    P3 = spdiagm([1.0, 0.0, 0.0])
    b3 = 0.0 #-4.0
    constr3 = nl_quad(a3,P3,b3)

    # x^2 + y^2 - 6x - 6y - s + 17 = 0
    # (x-3)^2 + (y-3)^2 - s - 1 = 0
    # (x-3)^2 + (y-3)^2 <= 1
    # ball function?
    a3 = sparse([-6.0, -6.0, -1.0])
    P3 = spdiagm([1.0, 1.0, 0.0])
    b3 = -17.0 #-4.0
    constr4 = nl_quad(a3,P3,b3)

    # convex problem
    println("CONVEX INFEASIBLE PROBLEM")
    nlp = class_nlp(3, obj, [constr2, constr4])
    vars = class_variables(3,2);
    status, vars = homogeneous_algorithm(nlp, vars, settings)
    @show x_scaled(vars)
    @show y_scaled(vars)

    println("NONCONVEX INFEASIBLE PROBLEM 1")
    nlp = class_nlp(3, obj, [constr1, constr2])
    vars = class_variables(3,2);
    status, vars = homogeneous_algorithm(nlp, vars, settings)
    @show x_scaled(vars)
    @show y_scaled(vars)

    println("NONCONVEX INFEASIBLE PROBLEM 2")
    nlp = class_nlp(3, obj, [constr1_b, constr2])
    vars = class_variables(3,2);
    status, vars = homogeneous_algorithm(nlp, vars, settings)
    @show x_scaled(vars)
    @show y_scaled(vars)

    nlp = class_nlp(3, obj, [constr1, constr3])
    vars = class_variables(3,2);
    status, vars = homogeneous_algorithm(nlp, vars, settings)
    @show x_scaled(vars)
    @show y_scaled(vars)
    #@show internal_eval_hesslag_prod(nlp,x_scaled(vars),y_scaled(vars))

    nlp = class_nlp(3, obj, [constr2, constr3])
    vars = class_variables(3,2);
    status, vars = homogeneous_algorithm(nlp, vars, settings)
    @show x_scaled(vars)
    @show y_scaled(vars)

    #@fact status --> :locally_optimal
    #@fact norm(x_scaled(vars) - [2.0,2.0,0.0]) --> less_than(1e-6)
    #@fact tau(vars) --> greater_than(0.1)
    #@fact 0.0 --> less_than(kappa(vars))
    #@fact kappa(vars) --> less_than(1e-6)
end
