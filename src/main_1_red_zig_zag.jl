using MeshCat
using Rotations
using ForwardDiff
using ControlSystems

include("visualize.jl")
include("quaternion.jl")
include("quadrotor_dynamics.jl")
include("lqr.jl")

blue = Quadrotor()
red = Quadrotor()

times = LinRange(0, 10, 501)
dt = times[2]-times[1]

red_cost_Q̃ = Diagonal([10, 10, 10, 1, 1, 1, 10, 10, 10, 1, 1, 1])
red_cost_R = Diagonal([0.1, 0.1, 0.1, 0.1])
uhover = [0.0, 0.0, 0.0, 0.0] .+ red.m*red.g/4.0

x0_blue = [0.0, 5.0, 14.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
Xref_blue = [x0_blue for i=1:length(times)]

u0_red = copy(uhover)
Uref_red = [copy(u0_red) for i=1:length(times)-1]

x0_red = [ 0.0,  100.0,  100.0,  1.0,  0.0,  0.0,  0.0,  0.0,  -10,  -10,  0.0,  0.0,  0.0]
Xref_red = [copy(x0_red) for i=1:length(times)]
for i=1:length(times)-1
  Xref_red[i][10] = copy(Xref_red[i][10]) + 20*cos(pi/32*i)
  xref = copy(Xref_red[i])
  Xref_red[i+1][1:3] = xref[1:3] + Xref_red[i][8:10] * dt
end
# for i=126:250
#   xref = copy(Xref_red[i])
#   Xref_red[i][8:10] = [0.0, 0.0, -20.0]
#   Xref_red[i+1][1:3] = xref[1:3] + Xref_red[i][8:10]*dt
# end
# for i=251:375
#   xref = copy(Xref_red[i])
#   Xref_red[i][8:10] = [0.0, -20.0, 0.0]
#   Xref_red[i+1][1:3] = xref[1:3] + Xref_red[i][8:10]*dt
# end
# for i=376:500
#   xref = copy(Xref_red[i])
#   Xref_red[i][8:10] = [0.0, 0.0, -20.0]
#   Xref_red[i+1][1:3] = xref[1:3] + Xref_red[i][8:10]*dt
# end

A = [zeros(12,12) for i=1:length(times)-1]
B = [zeros(12,4) for i=1:length(times)-1]

for i=1:length(times)-1
  A[i], B[i] = get_Ã_B̃(red, Xref_red[i], Xref_red[i+1], Uref_red[i], red_cost_Q̃, red_cost_R, dt)
end
red_K, = riccati(A, B, red_cost_Q̃, red_cost_R, red_cost_Q̃)

red_lqr = LQRController(red_K, times, Xref_red, Uref_red)

X_red = [copy(x0_red) for i=1:length(times)]
U_red = [copy(uhover) for i=1:length(times)-1]

for i=1:length(times)-1
  U_red[i] = get_control(red_lqr, X_red[i], i)
  X_red[i+1] = quad_dynamics_rk4(red, X_red[i], U_red[i], times[i+1]-times[i])
end

# vis = initialize_visualizer()
visualize!(vis, times[end], Xref_blue, X_red)
visualize!(vis, times[end], Xref_blue, Xref_red)


#############
x_red_next = copy(X_red[1])
test_i=1

u_red_next = copy(get_control(red_lqr, x_red_next, test_i))
x_red_next = copy(quad_dynamics_rk4(red, x_red_next, u_red_next, dt))
test_i=test_i+1
red_K[4]