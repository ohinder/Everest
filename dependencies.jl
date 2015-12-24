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

using MPI
x=1+1

import MUMPS # needs to be intialized before Ipopt don't understand why
MPI.Init();
cntl = MUMPS.default_icntl;
cntl[4] = 1;
var_mumps = MUMPS.Mumps{Float64}(MUMPS.mumps_unsymmetric, cntl, MUMPS.default_cntl64	);  # Real, general unsymmetric
A = sparse(rand(4,4)); rhs = rand(4);       # Happens on all cores
MUMPS.associate_matrix(var_mumps, A);
MUMPS.factorize(var_mumps);
MUMPS.associate_rhs(var_mumps, rhs);
MUMPS.solve(var_mumps);
x_sol = MUMPS.get_solution(var_mumps);
println(var_mumps.infog[12])
MUMPS.finalize(var_mumps);
MPI.Finalize();

Pkg.update()



Pkg.clone("https://github.com/ohinder/Everest.git")
Pkg.clone("https://github.com/ohinder/my_advanced_timer.git")
Pkg.update()

using Everest
