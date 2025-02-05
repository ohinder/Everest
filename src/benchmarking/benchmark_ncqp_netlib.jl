export test_problem_ncqp, tridiagonal

#
# test the linear program solver on a simple problem
#
using Ipopt, KNITRO
#simple_tests();

function test_problem_ncqp(name::String)
	file_name = name * ".mat";

  dir = dirname(@__FILE__) * "/Problems/"

  A, b, c = get_netlib_problem(dir, file_name);
  Q = tridiagonal(length(c), 0.0,  1.0);

	println("Solving ", file_name, " with the homogeneous algorithm")
	println(size(A,2), " variables and ", size(A,1), " constraints")
	println("Non-zeros: ", length(nonzeros(A)))
	settings = class_settings();

	#println("=================== Linear system solver is matlab ldl ==================")
	println("=================== Linear system solver is MUMPS ==================")
	settings.linear_system_solver = linear_solver_MUMPS();#linear_solver_MATLAB();

	settings.newton_solver = class_homogeneous_newton();

  qp = class_quadratic_program(A, b, c, Q);
  homogeneous_algorithm(qp,settings)
  #trivial_test(A, b, c, Q, 1, file_name, settings, true);

	#settings.newton_solver = class_newton_ip();
	#trivial_test(A, b, c, Q, 1, file_name, settings, true);

	println("=========================================================================")
	println("Calls IPOPT")
	solve_with_JuMP(A, b, c, Q, IpoptSolver(max_iter=1000,print_level=3));
end
