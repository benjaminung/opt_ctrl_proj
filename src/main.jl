using MeshCat
using Rotations
using ForwardDiff
using ControlSystems

include("visualize.jl")
include("quaternion.jl")
include("quadrotor_dynamics.jl")
include("lqr.jl")
# test_x = zeros(13)
# test_x[4] = 1.0
# test_u = zeros(4)
# sum(RobotZoo.trim_controls(model))
# model_blue.mass * model_blue.gravity
blue = Quadrotor()
red = Quadrotor()

times = LinRange(0, 1, 51)
dt = times[2]-times[1]

red_cost_Q̃ = Diagonal([10, 10, 10, 1, 1, 1, 10, 10, 10, 1, 1, 1])
red_cost_R = Diagonal([0.1, 0.1, 0.1, 0.1])

x0_blue = [0.0, 5.0, 14.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
x0_red = [0.0, 100.0, 100.0, 1.0, 0.0, 0.0, 0.0, 0.0, -10.0, -10.0, 0.0, 0.0, 0.0]

Xref_blue = [copy(x0_blue) for i=1:length(times)]
Uref_blue = [zeros(4) for i=1:length(times)]

# xeq_red = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -10.0, -10.0, 0.0, 0.0, 0.0]
# A, B = get_Ã_B̃(red, xeq_red, ueq_red, red_cost_Q̃, red_cost_R, dt) 
ueq_red = [0.0, 0.0, 0.0, 0.0] #.+ red.m*red.g/4.0

Xref_red = [copy(x0_red) for i=1:length(times)]
Xeq_red = [copy(xeq_red) for i=1:length(times)]
Ueq_red = [copy(ueq_red) for i=1:length(times)]

xeq_1_red = [-20.0, 60.0, 100.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
A1, B1 = get_Ã_B̃(red, xeq_1_red, ueq_red, red_cost_Q̃, red_cost_R, dt)
# xeq_2_red = [-20.0, 60.0, 60.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -20.0, 0.0, 0.0, 0.0]
# A2, B2 = get_Ã_B̃(red, xeq_2_red, ueq_red,red_cost_Q̃, red_cost_R, dt)
# xeq_3_red = [0.0, 60.0, 60.0, 1.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# A3, B3 = get_Ã_B̃(red, xeq_3_red, ueq_red,red_cost_Q̃, red_cost_R, dt)

red_K = [zeros(size(B1')) for i=1:length(times)]

Uref_red = [zeros(4) for i=1:length(times)]
for i=1:50
  Xref_red[i][8:10] .= [-20.0, -40.0, 0.0]
  Xeq_red[i] = [-20.0, 60.0, 100.0, 1.0, 0.0, 0.0, 0.0, -20.0, -40.0, 0.0, 0.0, 0.0, 0.0]
end
red_K[1:50], = riccati(A1, B1, red_cost_Q̃, red_cost_R, red_cost_Q̃, 51)

for i=51:150
  Xref_red[i][8:10] .= [0.0, 0.0, -20.0]
  Xeq_red[i] = [-20.0, 60.0, 60.0, 1.0, 0.0, 0.0, 0.0, -20.0, 0.0, 0.0, 0.0, 0.0]
end
red_K[51:150], = riccati(A2, B2, red_cost_Q̃, red_cost_R, red_cost_Q̃, 101)

for i=151:200
  Xref_red[i][8:10] .= [20.0, 0.0, 0.0]
  Xeq_red[i] = [0.0, 60.0, 60.0, 1.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 0.0]
end
red_K[151:200], = riccati(A3, B3, red_cost_Q̃, red_cost_R, red_cost_Q̃, 51)
red_K[201:501], = riccati(A, B, red_cost_Q̃, red_cost_R, red_cost_Q̃, 302)

for i=1:length(times)-1
  Xref_red[i+1][1:3] .= Xref_red[i][1:3] + Xref_red[i][8:10] * (times[i+1]-times[i])
  Xref_blue[i+1] = quad_dynamics_rk4(blue, Xref_blue[i], [1.23625, 1.23625, 1.23625, 1.23625], times[i+1]-times[i])
end

red_lqr = LQRController(red_K, times, Xeq_red, Ueq_red, false)

X_red = copy(Xref_red)
U_red = copy(Uref_red)


x_red_next = copy(X_red[1])
test_i=1

u_red_next = copy(get_control(red_lqr, x_red_next, test_i))
x_red_next = copy(quad_dynamics_rk4(red, x_red_next, u_red_next, dt))
test_i=test_i+1
red_K[4]

for i=1:length(times)-1
  U_red[i] = get_control(red_lqr, X_red[i], i)
  X_red[i+1] = quad_dynamics_rk4(red, X_red[i], U_red[i], times[i+1]-times[i])
end
# Xref_blue = copy(Xref)

u0 = [0.0, 0.0, 0.0, 0.0]
Uref = [u0 for i=1:length(Xref_blue)-1]

#how to rotate
# let
#   test_rot_ang = pi/4.0
#   test_quat = [cos(test_rot_ang/2.0), sin(test_rot_ang/2.0), 0.0, 0.0]
#   test_pos = [0.0, 0.0, 1.0]
#   rot = UnitQuaternion(test_quat)
#   rot * test_pos
# end

# vis = Visualizer()
# set_mesh!(vis, model)
# render(vis)
# visualize!(vis, model, times[end], Xref)

vis = initialize_visualizer()
# render(vis)
visualize!(vis, times[end], Xref_blue, X_red)
visualize!(vis, times[end], Xref_blue, Xref_red)

# model_rocker = RobotZoo.PlanarRocket()
# rocke_x0 = zeros(8)
# rocke_u0 = zeros(2)
# size(model_rocker)
# discrete_dynamics(RK4, model_rocker, rocke_x0, rocke_u0, times[1], times[2]-times[1])

# RobotDynamics.DiscretizedDynamics(quadrotors.model_blue, RobotDynamics.RK4)

# x_1 = discrete_jacobian(RK4, quadrotors.model_blue, SVector{13,Float64}(Xref_blue[1]), SVector{4,Float64}([2, 2, 1, 1]), times[1], times[2]-times[1])
# zeq = KnotPoint(SVector{13,Float64}(Xref_blue[1]), SVector{4,Float64}(Uref[1]), times[2]-times[1])
# ∇f = RobotDynamics.DynamicsJacobian(quadrotors.model_blue)
# jacobian!(∇f, quadrotors.model_blue, zeq)
# x_1 = [Xrefxdot[1:3]; xdot[4:7]
# RobotDynamics.RK4(quadrotors.model_blue)
# RobotDynamics.QuadratureRule.RK4

cd("/OptCtrl/Project/optimal_control_project/src")
import Pkg; Pkg.activate(joinpath(@__DIR__,"..")); Pkg.instantiate()

model_e = Quadrotor{Rotation{3}}()
size(model_e)
size(quadrotors.model_blue)

using LinearAlgebra
using SparseArrays
using ForwardDiff
using OSQP
using RobotDynamics
using RobotZoo: Quadrotor
using RobotZoo
using StaticArrays
using Plots
using Statistics

using MeshCat
using CoordinateTransformations
using Rotations 
using Colors
using GeometryBasics
using RobotZoo: PlanarRocket

# model = PlanarRocket(max_roll=5.0)
# const n = state_dim(model)
# const m = control_dim(model)
# const xeq = [zeros(6); model.m*model.g; 0]
# const ueq = zeros(2) 
# norm(dynamics(model, xeq, ueq)) ≈ 0 # make sure it's an equilibrium point
# const dt = 0.1  # time step (s)
# const tf = 25   # time horizon (s)
# const N = Int(tf / dt) + 1

# zeq = KnotPoint(xeq,ueq,dt)   # create a `KnotPoint` type that stores everything together
# ∇f = RobotDynamics.DynamicsJacobian(model)
# jacobian!(∇f, model, zeq)
# discrete_jacobian!(RK4, ∇f, model, zeq)
# discrete_dynamics(RK4, model, SVector{8, Float64}([0, 0, 0, 0, 0, 0, 0, 0]), SVector{2, Float64}([0, 0]), times[1], times[2]-times[1])