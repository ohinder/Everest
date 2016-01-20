info("Loading required packages")
using advanced_timer, MATLAB, MUMPS

info("Loading julia code")

include("linear_system_solvers/linear_system_solvers.jl")
include("core/core.jl")
include("helpers/eigenvalues.jl")
include("helpers/woodbury.jl")
include("newton_direction/newton_direction.jl")
include("line_search/line_search.jl")
include("strategies/strategies.jl")
include("homogeneous.jl")
include("benchmarking/benchmarking.jl")
