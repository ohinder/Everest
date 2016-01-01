# About Everest

* work in progress *

Everest is an interior point optimization solver for problems of the form:

min f(x)

Ax = b

x >= 0

Where f(x) is any smooth non-convex function.

# How to install

You need to install at least one of the following package manually:

[MATLAB](https://github.com/JuliaLang/MATLAB.jl)

[MUMPS](https://github.com/JuliaOptimizers/MUMPS.jl)

Currently these packages are only supported in Linux or Mac OS X.

After installing these packages run:

```julia
Pkg.add("JuMP")
Pkg.add("MAT")
Pkg.clone("https://github.com/ohinder/Everest.jl.git")
```

To test your installation run:

```julia
using Everest
Pkg.test("Everest")
```
