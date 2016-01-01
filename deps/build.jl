run(`make MUMPS`)
include("deps.jl")

Pkg.build("MUMPS")
