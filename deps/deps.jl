Pkg.add("Ipopt")
#Pkg.add("KNITRO")
Pkg.add("JuMP")
Pkg.add("MAT")
Pkg.clone("https://github.com/JuliaOptimizers/MUMPS.jl")

Pkg.update()
