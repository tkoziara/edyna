using LinearAlgebra
include("../ellell.jl")
include("../expmap.jl")

x1 = [0.; 0.; 0.];
a1 = 1.;
b1 = 1.;
c1 = 1.;
R1 = Matrix{Float64}(I, 3, 3)

x2 = [0.; 0.; 2.0];
a2 = 1.;
b2 = 1.;
c2 = 1.1;
R2 = Matrix{Float64}(I, 3, 3)

(conpnt, normal, depth, A, B) = ellipsoid_ellipsoid_contact(x1, a1, b1, c1, R1, x2, a2, b2, c2, R2)

println("-------------------------------------------------------------")
println("for x1 = ", x1, "; a1, b1, c1 = ", a1, ", ", b1, ", ", c1, "; R1 = ", R1)
println("and x2 = ", x2, "; a2, b2, c2 = ", a2, ", ", b2, ", ", c2, "; R2 = ", R2)
println("conpnt = ", conpnt, ", normal = ", normal, ", depth = ", depth, ", A, B = ", A, ", ", B)
