cd("src")
include("visualize.jl")
include("quaternion.jl")
include("quadrotor_dynamics.jl")
include("lqr.jl")
include("mpc.jl")
include("red_traj.jl")

blue = Quadrotor()
red = Quadrotor()

times = LinRange(0, 20, 1001)

blue_cost_Q̃ = Diagonal([10, 10, 10, 1, 1, 1, 10, 10, 10, 1, 1, 1])
blue_cost_R = Diagonal([0.1, 0.1, 0.1, 0.1])

x0_blue = [0.0, -5.0, 14.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
xeq_blue = [0.0, 14.142, 14.142, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
Xref_blue = [copy(x0_blue) for i=1:length(times)]

uhover_blue = [0.0, 0.0, 0.0, 0.0] .+ blue.m*blue.g/4.0
u0_blue = copy(uhover_blue)
Uref_blue = [copy(u0_blue) for i=1:length(times)-1]

mpc_blue = OSQPController(n, m, Nmpc, length(Xref))
mpc_blue.times .= times
mpc_blue.Xref .= Xref_blue
mpc_blue.Uref .= Uref_blue

update_xref!(mpc_blue, x0, times[1], dt, xeq, times[end])

A_blue, B_blue = get_Ã_B̃(red, mpc_blue.Xref[end], mpc_blue.Xref[end], mpc_blue.Uref[end], blue_cost_Q̃, blue_cost_R, dt)
buildQP_constrained!(mpc_blue, A_blue, B_blue, blue_cost_Q̃, blue_cost_R, blue_cost_Q̃)


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