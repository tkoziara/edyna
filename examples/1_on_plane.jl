# 1 ellipsoid on a plane

#push!(LOAD_PATH, "/Users/tomek/Dropbox/tkwork/codes/edyna")
#import edyna
include("../edyna.jl")

edyna.ellipsoid([0, 0, 0.1], 0.01, 0.02, 0.03, [1. 0. 0.; 0. 1. 0.; 0. 0. 1.], 1000)

edyna.obstacle([[[-1, -1, 0], [-1, 1, 0], [1, -1, 0]], [[1, -1, 0], [-1, 1, 0], [1, 1, 0]]])

edyna.gravity(0, 0, -10)

edyna.run_rigid(1.0, 0.001, 0.01, "out/1_on_plane")
