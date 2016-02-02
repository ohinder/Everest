using MPI
import MUMPS # needs to be intialized before Ipopt don't understand why

type linear_solver_MUMPS <: abstract_linear_system_solver
	_SparseMatrix::SparseMatrixCSC{Float64,Int64}
	_factor::MUMPS.Mumps{Float64} # TO DO, give type

  function linear_solver_MUMPS()
      return new();
  end
end

# intialize and finalize function???
# MPI ???
function mumps_sym(sym::Symbol)
	if sym == :symmetric
			return MUMPS.mumps_symmetric
	elseif sym == :unsymmetric
			return MUMPS.mumps_unsymmetric
	elseif sym == :definite
			return MUMPS.mumps_definite
	else
			error("this symmetry symbol is not understood")
	end
end

function initialize!(solver::linear_solver_MUMPS, A::SparseMatrixCSC{Float64,Int64}, sym::Symbol)
		if ~MPI.Initialized()
			MPI.Init()
		end

		icntl = MUMPS.default_icntl[:]; # copy
    icntl[4] = 1;
		icntl[14] = 100.0;

    solver._factor = MUMPS.Mumps{Float64}(mumps_sym(sym), icntl, MUMPS.default_cntl64	);  # Real, general unsymmetric
		solver._SparseMatrix = A;

    return solver;
end

function finalize!(solver::linear_solver_MUMPS)
    MUMPS.finalize(solver._factor);
end

function ls_factor!(solver::linear_solver_MUMPS, n::Int64, m::Int64)
		start_advanced_timer("MUMPS/associate_matrix")
    MUMPS.associate_matrix(solver._factor, solver._SparseMatrix);
		pause_advanced_timer("MUMPS/associate_matrix")
		
		MUMPS.factorize(solver._factor);

		return solver._factor.infog[12] == m; # inertia
end

function ls_solve(solver::linear_solver_MUMPS, my_rhs::AbstractArray)
    MUMPS.associate_rhs(solver._factor, my_rhs);
    MUMPS.solve(solver._factor);
    return MUMPS.get_solution(solver._factor);
end
