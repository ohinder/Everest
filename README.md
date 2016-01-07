# About Everest

**This project is work in progress.**

Everest is an interior point optimization solver for problems of the form:

min f(x)

Ax = b

x >= 0

Where f(x) is any smooth non-convex function.

# How to install

1. Install at least one of the following package manually (the linear system solver):

[MATLAB.jl](https://github.com/JuliaLang/MATLAB.jl)

[MUMPS.jl](https://github.com/JuliaOptimizers/MUMPS.jl)

MUMPS is preferable, since it is open source and is significantly faster than using the MATLAB LBL factorization (less overhead). If you want to add an new linear solver you can do so in this [directory](https://github.com/ohinder/Everest.jl/tree/master/src/linear_system_solvers).

Currently these packages are only supported in Linux or Mac OS X.

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
