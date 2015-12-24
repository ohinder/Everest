using Base.Test
# test linear system solvers

function test_linear_solver(solver::abstract_linear_system_solver, sym::Symbol)
    n = 40;
    mat = speye(n)
    mat = mat + mat'
    mat[1,1] = -100;
    mat[2,2] = -100;

    initialize!(solver, mat, sym);

    test_linear_solver_inertia(solver, sym, mat, 1, 38, 2)

    rhs = 1.0*rand(n, 2);
    A = rand(n,n)
    mat[:,:] = A + A'

    test_linear_solver_accuracy(solver, sym, mat, rhs, 38, 2)

    finalize!(solver)
end

function test_linear_solver_inertia(solver::abstract_linear_system_solver, sym::Symbol, mat::SparseMatrixCSC{Float64,Int64}, correct_inertia::Int64, n::Int64, m::Int64)

   @test solver._SparseMatrix == mat

   inertia = ls_factor!(solver, n, m)
   @test inertia == correct_inertia
end

function test_linear_solver_accuracy(solver::abstract_linear_system_solver, sym::Symbol, mat::SparseMatrixCSC{Float64,Int64}, rhs::AbstractArray, n::Int64, m::Int64)
  ls_factor!(solver, n, m)

  sol = ls_solve(solver,rhs);

  #sol = ls_solve(solver,rhs);
  err = norm(mat*sol - rhs,1)/norm(sol,1)
  @test err < 1e-6
  @show err
end

# test mumps solver
begin
  ls_solver_mumps = linear_solver_MUMPS();
  test_linear_solver( ls_solver_mumps, :symmetric )
  #MPI.Finalize()
end

# test julia solver
begin
  # test julia lu factor
  ls_solver_julia = linear_solver_JULIA();
  test_linear_solver( ls_solver_julia , :unsymmetric )

  # test julia cholesky factor
  ls_solver_julia = linear_solver_JULIA();
  test_linear_solver( ls_solver_julia, :definite )
end

# test matlab ldl solver
begin
  ls_solver_matlab = linear_solver_MATLAB();
  test_linear_solver( ls_solver_matlab, :symmetric )
end
