Pkg.add("Ipopt")
Pkg.add("KNITRO")
Pkg.add("JuMP")

Pkg.clone("https://github.com/JuliaLang/MATLAB.jl.git")

# install http://brew.sh
# brew tap homebrew/science
# brew install mumps [--with-scotch5] [--with-openblas]  # use scotch/openblas at your option
Pkg.clone("https://github.com/JuliaOptimizers/MUMPS.jl.git")
Pkg.build("MUMPS")

using MPI, MUMPS
let
  MPI.Init();
  #mumps = Mumps{Float64}(mumps_definite);

  A = sprand(10, 10, .2) + speye(10); rhs = rand(10);
  x = solve(A, rhs);  # Mumps object is created and destroyed
  norm(x - A \ rhs) / norm(x)
  MPI.Finalize();
end

MPI.Init();
mumps = Mumps{Float64}(mumps_unsymmetric, default_icntl, default_cntl64	);  # Real, general unsymmetric
A = sparse(rand(4,4)); rhs = rand(4);       # Happens on all cores
associate_matrix(mumps, A);
factorize(mumps);
associate_rhs(mumps, rhs);
solve(mumps);
x = get_solution(mumps);
finalize(mumps);
MPI.Finalize();

Pkg.update()



Pkg.clone("https://github.com/ohinder/Everest.git")
Pkg.clone("https://github.com/ohinder/my_advanced_timer.git")
Pkg.update()

using Everest
