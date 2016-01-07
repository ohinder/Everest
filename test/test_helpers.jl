function test_woodbury()
    len = 500;
    ls_solver_julia = linear_solver_JULIA();
    solver = ls_solver_julia;
    mat = sparse(rand(len,len));
    initialize!(solver, mat, :unsymmetric);

    @test solver._SparseMatrix == mat

    inertia = ls_factor!(solver, len, 0)
    @test inertia == 1

    sol = 1.0*zeros(len);
    rhs = 1.0*ones(len);

    U = rand(len,2);
    #U[1,1] = 2.0;
    #U[2,2] = 2.0;
    V = copy(U)';

    W = woodbury_identity(solver, U, V);

    sol = ls_solve(W,rhs);

    err = norm((mat + U * V) * sol - rhs)
    @show err
    @fact err --> less_than(1e-10)
end

facts("helpers") do
  context("woodbury") do
      test_woodbury();
  end
end
