# Two levels of tests
# Tests of individual functions
# Tests of code to solve simple unbounded, bounded problems etc ...

using Base.Test
using FactCheck

begin
  # includ
  include("../src/loadcode.jl");

  include("test_linear_system_solvers.jl")
  include("test_helpers.jl")
  include("test_algorithm_components.jl") # merge tests into test algorithm components
  include("sample_problem_tests.jl")
end
