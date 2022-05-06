using MeshCat
using Rotations
using ForwardDiff
using ControlSystems

include("visualize.jl")
include("quaternion.jl")
include("quadrotor_dynamics.jl")
include("lqr.jl")
include("red_traj.jl")

blue = Quadrotor()
red = Quadrotor()

times = LinRange(0, 20, 1001)

x0_blue = [0.0, 5.0, 14.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
Xref_blue = [x0_blue for i=1:length(times)]

X_red = red_traj(red, times, zig_zag_z, 0.2, 6)

vis = initialize_visualizer()
visualize!(vis, times[end], Xref_blue, X_red)

#############
x_red_next = copy(X_red[1])
test_i=1

u_red_next = copy(get_control(red_lqr, x_red_next, test_i))
x_red_next = copy(quad_dynamics_rk4(red, x_red_next, u_red_next, dt))
test_i=test_i+1
red_K[4]