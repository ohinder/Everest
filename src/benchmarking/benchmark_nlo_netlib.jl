export test_problem, tridiagonal

#
# test the linear program solver on a simple problem
#
using Ipopt, KNITRO
#simple_tests();

# 54
# 59
#

function solve_with_JuMP_sqrt(A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, c::Array{Float64,1}, solver=IpoptSolver(max_iter=300))
	model = Model(solver=solver)

	n, k = size(A)

	@defVar(model, x[1:k] >= 0)
	@setNLObjective( model, Min, sum{log(x[i]+1), i=1:k})


	#@addConstraint(model, A*x .== b)
	@defConstrRef constr[1:n]
	for i = 1:n
		constr[i] = @addConstraint(model, sum{A[i,j]*x[j], j=1:k} == b[i])
	end

	status = solve(model)

	if status == :Optimal
		return 1, getValue(x)[:], getDual(constr)[:]
	elseif status == :Infeasible
		return 2
	elseif status == :Unbounded
		return 3
	else
		return 0
	end
end

function test_problem_nlo(name::String)
	file_name = name * ".mat";

  dir = dirname(@__FILE__) * "/Problems/"

  A, b, c = get_netlib_problem(dir, file_name);
  Q = tridiagonal(length(c),0.0, 0.0);
  #Q = tridiagonal(length(c),-1.0,0.0);

	println("Solving ", file_name, " with the homogeneous algorithm")
	println(size(A,2), " variables and ", size(A,1), " constraints")
	println("Non-zeros: ", length(nonzeros(A)))
	settings = class_settings();

	println("=================== Linear system solver is MUMPS ==================")
	settings.linear_system_solver = linear_solver_MUMPS();

	settings.newton_solver = class_homogeneous_newton();

  qp = class_quadratic_program(A, b, c, Q);
  nlo = class_nlo(A, b, nl_log(length(c)))
  homogeneous_algorithm(nlo,settings)
  #trivial_test(A, b, c, Q, 1, file_name, settings, true);

	#settings.newton_solver = class_newton_ip();
	#trivial_test(A, b, c, Q, 1, file_name, settings, true);

	println("=========================================================================")
	println("Calls IPOPT")
	solve_with_JuMP_sqrt(A, b, c, IpoptSolver(max_iter=1000,print_level=3));
end
