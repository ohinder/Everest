# About Everest

**This project is work in progress.**

Everest is an interior point optimization solver for problems of the form:

min f(x)

Ax = b

x >= 0

Where f(x) is any smooth non-convex function.

# How to install

1. Install the linear solver [MUMPS.jl](https://github.com/JuliaOptimizers/MUMPS.jl). Currently this package is only supported in Linux or Mac OS X. If you want to add an new linear solver you can do so in this [directory](https://github.com/ohinder/Everest.jl/tree/master/src/linear_system_solvers).


2. Run the following code in julia:

```julia
Pkg.add("JuMP")
Pkg.add("MAT")
Pkg.clone("https://github.com/ohinder/advanced_timer.jl.git")
Pkg.clone("https://github.com/ohinder/Everest.jl.git")
```

3. To test your installation run in julia:

```julia
using Everest
Pkg.test("Everest")
```

# How to use
