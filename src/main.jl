import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate();
cd("src")
include("visualize.jl")
include("quaternion.jl")
include("quadrotor_dynamics.jl")
include("lqr.jl")
include("mpc.jl")
include("red_traj.jl")

blue_quad = Quadrotor()
red_quad = Quadrotor()

times = LinRange(0, 20, 1001)

X_red_const_vel = red_traj(red_quad, times, const_vel)
X_red_zigz = red_traj(red_quad, times, zig_zag_z, 0.2, 6)
X_red_zigy = red_traj(red_quad, times, zig_zag_y, 0.2, 5)
X_red_zigx = red_traj(red_quad, times, zig_zag_x, 0.2, 5)

vis = initialize_visualizer()
# visualize!(vis, times[end], Xref_blue, X_red) # visualize X_red

x0_blue = [0.0, -5.0, 14.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
xeq_blue = [0.0, 20.0, 20.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
Xref_blue = [copy(xeq_blue) for i=1:length(times)]

uhover_blue = [0.0, 0.0, 0.0, 0.0] .+ blue_quad.m*blue_quad.g/4.0
u0_blue = copy(uhover_blue)
Uref_blue = [copy(u0_blue) for i=1:length(times)]

blue_cost_Q̃ = Diagonal([10, 10, 10, 1, 1, 1, 10, 10, 10, 1, 1, 1])
blue_cost_R = Diagonal([0.1, 0.1, 0.1, 0.1])

n = 12 # number of states in mpc x̃
m = 4 # nummer of controls
dt = times[2]-times[1]
mpc_N = 20 # mpc horizon length
mpc_Nd = (mpc_N - 1)*(n + 8) # mpc constraints
mpc_blue = OSQPController(n, m, mpc_N, length(Xref_blue), mpc_Nd)
mpc_blue.times .= times
mpc_blue.Xref .= Xref_blue
mpc_blue.Uref .= Uref_blue

mpc_blue_const_vel = deepcopy(mpc_blue)
mpc_blue_zigz = deepcopy(mpc_blue)
mpc_blue_zigy = deepcopy(mpc_blue)
mpc_blue_zigx = deepcopy(mpc_blue)

X_blue_const_vel, U_blue_const_vel, times_const_vel, X_red_sim_const_vel = simulate(
  blue_quad, x0_blue, mpc_blue_const_vel, blue_cost_Q̃, blue_cost_R, blue_cost_Q̃, X_red_const_vel
);
visualize!(vis, times_const_vel[end], X_blue_const_vel, X_red_sim_const_vel)

X_blue_zigz, U_blue_zigz, times_zigz, X_red_sim_zigz = simulate(
  blue_quad, x0_blue, mpc_blue_zigz, blue_cost_Q̃, blue_cost_R, blue_cost_Q̃, X_red_zigz
);
visualize!(vis, times_zigz[end], X_blue_zigz, X_red_sim_zigz)

X_blue_zigy, U_blue_zigy, times_zigy, X_red_sim_zigy = simulate(
  blue_quad, x0_blue, mpc_blue_zigy, blue_cost_Q̃, blue_cost_R, blue_cost_Q̃, X_red_zigy
);
visualize!(vis, times_zigy[end], X_blue_zigy, X_red_sim_zigy)

X_blue_zigx, U_blue_zigx, times_zigx, X_red_sim_zigx = simulate(
  blue_quad, x0_blue, mpc_blue_zigx, blue_cost_Q̃, blue_cost_R, blue_cost_Q̃, X_red_zigx
);
visualize!(vis, times_zigx[end], X_blue_zigx, X_red_sim_zigx)
