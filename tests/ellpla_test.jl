using LinearAlgebra
include("../ellpla.jl")
include("../expmap.jl")

x = [0.; 0.; 1.];
a = 1.;
b = 1.;
c = 1.05;
R = Matrix{Float64}(I, 3, 3)
point = [0.; 0.; 0.];
normal = [0.; 0.; 1.];

(conpnt, depth, A, B) = ellipsoid_plane_contact(x, a, b, c, R, point, normal)

println("-------------------------------------------------------------")
println("for x = ", x, "; a, b, c = ", a, ", ", b, ", ", c, "; R = ", R)
println("and point = ", point, ", normal = ", normal)
println("conpnt = ", conpnt, ", depth = ", depth, ", A, B = ", A, ", ", B)

a = 0.1
x[3] = 1.049

(conpnt, depth, A, B) = ellipsoid_plane_contact(x, a, b, c, R, point, normal)

println("-------------------------------------------------------------")
println("for x = ", x, "; a, b, c = ", a, ", ", b, ", ", c, "; R = ", R)
println("and point = ", point, ", normal = ", normal)
println("conpnt = ", conpnt, ", depth = ", depth, ", A, B = ", A, ", ", B)

b = 0.1
x[3] = 1.0499

(conpnt, depth, A, B) = ellipsoid_plane_contact(x, a, b, c, R, point, normal)

println("-------------------------------------------------------------")
println("for x = ", x, "; a, b, c = ", a, ", ", b, ", ", c, "; R = ", R)
println("and point = ", point, ", normal = ", normal)
println("conpnt = ", conpnt, ", depth = ", depth, ", A, B = ", A, ", ", B)

a = 0.5
b = 0.7
x[3] = 0.7
R = expmap([1.0 1.0 1.0])

(conpnt, depth, A, B) = ellipsoid_plane_contact(x, a, b, c, R, point, normal)

println("-------------------------------------------------------------")
println("for x = ", x, "; a, b, c = ", a, ", ", b, ", ", c, "; R = ", R)
println("and point = ", point, ", normal = ", normal)
println("conpnt = ", conpnt, ", depth = ", depth, ", A, B = ", A, ", ", B)
