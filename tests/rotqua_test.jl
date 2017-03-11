include("../expmap.jl")
include("../rotqua.jl")

R = eye(3);
q = rotqua(R)
println("-------------------------------------------------------------")
println("for R = ", R)
println("q = ", q)

R = expmap([1.0 1.0 1.0])
q = rotqua(R)
println("-------------------------------------------------------------")
println("for R = ", R)
println("q = ", q)


