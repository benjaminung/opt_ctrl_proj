include("visualize.jl")
include("point_dynamics.jl")
include("mpc.jl")

point = Point()

times_point = LinRange(0, 10, 501)
dt_point = times_point[2]-times_point[1]

point_cost_Q = Diagonal([10, 10, 10, 1, 1, 1])
point_cost_R = Diagonal([0.1, 0.1, 0.1])

x0_point = [
  -10.0, -10.0, 10.0, 
  0.0, 0.0, 0.0
  ]
u0_point = [0.0, 0.0, 9.81*point.m]
xeq_point = [
  10.0, 20.0, 20.0, 
  0.0, 0.0, 0.0
]

A = ForwardDiff.jacobian(_x->point_dynamics_rk4(point, _x, u0_point, dt_point), xeq_point)
B = ForwardDiff.jacobian(_u->point_dynamics_rk4(point, xeq_point, _u, dt_point), u0_point)

Xref_point = [copy(xeq_point) for i=1:length(times_point)]
Uref_point = [copy(u0_point) for i=1:length(times_point)]

n_point = 6
m_point = 3
mpc_point_N = 20
mpc_point_Nd = (mpc_point_N - 1)*(n_point+4)
mpc_point = OSQPController(n_point, m_point, mpc_point_N, length(Xref_point), mpc_point_Nd)
mpc_point.times .= times_point
mpc_point.Xref .= Xref_point
mpc_point.Uref .= Uref_point

buildQP_constrained_point!(point, mpc_point, A, B, point_cost_Q, point_cost_R, point_cost_Q)

X_point, U_point, = simulate_point(point, x0_point, mpc_point, point_cost_Q, point_cost_R, point_cost_Q, A)

X_point
U_point

x_point_next, u_point_calc, = simulate_one_step_point(point, x0_point, mpc_point, point_cost_Q, point_cost_R, point_cost_Q, A, 1)
x_point_next, u_point_calc, = simulate_one_step_point(point, x_point_next, mpc_point, point_cost_Q, point_cost_R, point_cost_Q, A, 2)