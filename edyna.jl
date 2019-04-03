include("ellpro.jl")

# edyna Julia module begins -->
module edyna

# ellipsoids data
ellcount = 0
ellcenter = [Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount)]
ellradii = [Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount)]
ellorient = [Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount),
             Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount),
             Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount)]
elldefo = [Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount),
           Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount),
           Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount)]
ellmass = Array{Float32}(undef, ellcount);
elleuler = [Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount),
            Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount),
            Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount)]
ellinveuler = [Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount),
               Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount),
               Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount)]
ellinertia = [Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount),
              Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount),
              Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount)]
ellinvinertia = [Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount),
                 Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount),
                 Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount)]
elllinear = [Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount)]
ellangular = [Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount)]

# insert ellipsoid
function ellipsoid(x, a, b, c, R, density, linear = [0.0 0.0 0.0], angular = [0.0 0.0 0.0])
  push!(ellcenter[1],x[1])
  push!(ellcenter[2],x[2])
  push!(ellcenter[3],x[3])
  push!(ellradii[1],a)
  push!(ellradii[2],b)
  push!(ellradii[3],c)
  push!(ellorient[1],R[1,1])
  push!(ellorient[2],R[2,1])
  push!(ellorient[3],R[3,1])
  push!(ellorient[4],R[1,2])
  push!(ellorient[5],R[2,2])
  push!(ellorient[6],R[3,2])
  push!(ellorient[7],R[1,3])
  push!(ellorient[8],R[2,3])
  push!(ellorient[9],R[3,3])
  push!(elldefo[1],1.0)
  push!(elldefo[2],0.0)
  push!(elldefo[3],0.0)
  push!(elldefo[4],0.0)
  push!(elldefo[5],1.0)
  push!(elldefo[6],0.0)
  push!(elldefo[7],0.0)
  push!(elldefo[8],0.0)
  push!(elldefo[9],1.0)

  (mass, volume, E, J) = ellipsoid_properties(a, b, c, R, density)

  push!(ellmass, mass)
  push!(elleuler[1], E[1,1])
  push!(elleuler[2], E[2,1])
  push!(elleuler[3], E[3,1])
  push!(elleuler[4], E[1,2])
  push!(elleuler[5], E[2,2])
  push!(elleuler[6], E[3,2])
  push!(elleuler[7], E[1,3])
  push!(elleuler[8], E[2,3])
  push!(elleuler[9], E[3,3])
  invE = inv(E)
  push!(ellinveuler[1], invE[1,1])
  push!(ellinveuler[2], invE[2,1])
  push!(ellinveuler[3], invE[3,1])
  push!(ellinveuler[4], invE[1,2])
  push!(ellinveuler[5], invE[2,2])
  push!(ellinveuler[6], invE[3,2])
  push!(ellinveuler[7], invE[1,3])
  push!(ellinveuler[8], invE[2,3])
  push!(ellinveuler[9], invE[3,3])
  push!(ellinertia[1], J[1,1])
  push!(ellinertia[2], J[2,1])
  push!(ellinertia[3], J[3,1])
  push!(ellinertia[4], J[1,2])
  push!(ellinertia[5], J[2,2])
  push!(ellinertia[6], J[3,2])
  push!(ellinertia[7], J[1,3])
  push!(ellinertia[8], J[2,3])
  push!(ellinertia[9], J[3,3])
  invJ = inv(J)
  push!(ellinvinertia[1], invJ[1,1])
  push!(ellinvinertia[2], invJ[2,1])
  push!(ellinvinertia[3], invJ[3,1])
  push!(ellinvinertia[4], invJ[1,2])
  push!(ellinvinertia[5], invJ[2,2])
  push!(ellinvinertia[6], invJ[3,2])
  push!(ellinvinertia[7], invJ[1,3])
  push!(ellinvinertia[8], invJ[2,3])
  push!(ellinvinertia[9], invJ[3,3])

  push!(elllinear[1],linear[1])
  push!(elllinear[2],linear[2])
  push!(elllinear[3],linear[3])
  push!(ellangular[1],angular[1])
  push!(ellangular[2],angular[2])
  push!(ellangular[3],angular[3])

  ellcount += 1
end

# triangles data
tricount = 0
trip = [Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount)]
triq = [Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount)]
trir = [Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount), Array{Float32}(undef, ellcount)]
trilinear = Array{Function}(undef, ellcount) # [vx, vy, vz] = linear[i](t)
triangular = Array{Function}(undef, ellcount) # [ox, oy, oz] = angular[i](t)

# insert an obstacle
zerofunc(t) = Float32[0.0; 0.0; 0.0]
function obstacle(p, q, r, linear = zerofunc, angular = zerofunc)
  for i = 1:size(p)[1]
    push!(trip,p[i])
    push!(triq,q[i])
    push!(trir,r[i])
    push!(trilinear, linear)
    push!(triangular, angular)
    tricount += 1
  end
end

# set gravity
grav = [0, 0, 0]
function gravity(g1, g2, g3)
 grav[1] = g1
 grav[2] = g2
 grav[3] = g3
end

end
# <-- edyna Julia module ends
