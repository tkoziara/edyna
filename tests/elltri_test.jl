using LinearAlgebra
include("../elltri.jl")
include("../expmap.jl")

x = [0.; 0.; 1.];
a = 1.;
b = 1.;
c = 1.05;
R = Matrix{Float64}(I, 3, 3)
p = [-1.; -1.; 0.];
q = [1.; -1.; 0.];
r = [0.; 1.; 0.];

(conpnt, normal, depth, A, B) = ellipsoid_triangle_contact(x, a, b, c, R, p, q, r)

println("-------------------------------------------------------------")
println("for x = ", x, "; a, b, c = ", a, ", ", b, ", ", c, "; R = ", R)
println("and p , q, r = ", p, ", ", q, ", ", r)
println("conpnt = ", conpnt, ", normal = ", normal, ", depth = ", depth, ", A, B = ", A, ", ", B)
