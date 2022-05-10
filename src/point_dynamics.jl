struct Point
  m::Float64
  g::Float64
  max_thrust_z::Float64
  min_thrust_z::Float64
  max_thrust_x_y::Float64
end

function Point()
  m = 0.5
  g = 9.81
  max_thrust_z = m*g*2.0
  min_thrust_z = 0
  max_thrust_x_y = m*10.0
  Point(m, g, max_thrust_z, min_thrust_z, max_thrust_x_y)
end

function point_dynamics(point::Point,x,u)
  accel = [0, 0, -point.g] + u/point.m
  vel = copy(x[4:6])
  return [vel; accel]
end

function point_dynamics_rk4(point::Point, x, u, h)
  f1 = point_dynamics(point, x, u)
  f2 = point_dynamics(point, x + 0.5*h*f1, u)
  f3 = point_dynamics(point, x + 0.5*h*f2, u)
  f4 = point_dynamics(point, x + h*f3, u)
  xn = Vector(copy(x) + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4))
  return xn
end

# function update_xref_point!(ctrl::MPCController, x, time, dt, xeq, time_eq)
#   t = get_k(ctrl, time)
#   # println(t)
#   teq = get_k(ctrl, time_eq)
#   # println(teq)
#   vel = (copy(xeq[1:3]) - copy(x[1:3]))/(time_eq-time)
#   # println(vel)
#   for k=t:teq - 1
#     ctrl.Xref[k][4:6] = copy(vel)
#     ctrl.Xref[k+1][1:3] = copy(ctrl.Xref[k][1:3]) + copy(ctrl.Xref[k][4:6])*dt
#   end
#   if teq < length(ctrl.times)
#     ctrl.Xref[teq][4:6] = [0.0, 0.0, 0.0]
#     for k=teq:length(ctrl.times)-1
#       ctrl.Xref[k+1][1:3] = copy(xeq[1:3])
#       ctrl.Xref[k+1][4:6] = [0.0, 0.0, 0.0]
#     end
#   end
# end

function get_control_point(ctrl::MPCController{OSQP.Model}, x, time, Q, R, Qf, A)
  # Update the QP
  updateQP_constrained_point!(ctrl, x, time, Q, R, Qf, A)
  OSQP.update!(ctrl.solver, q=ctrl.q, l=ctrl.lb, u=ctrl.ub)

  # Solve QP
  results = OSQP.solve!(ctrl.solver)
  Δu = results.x[1:3]
  println(Δu)
  
  k = get_k(ctrl, time)
  return ctrl.Uref[end] + Δu 
end

function buildQP_constrained_point!(point::Point, ctrl::MPCController, A,B,Q,R,Qf; kwargs...)
  # TODO: Implement this method to build the QP matrices
  n = length(ctrl.Xref[1])
  # SOLUTION:
  Nt = ctrl.Nmpc-1
  Nx = length(ctrl.Xref[1])    # number of states in x̃
  Nu = length(ctrl.Uref[1])    # number of controls

  # local ABs
  # for i=1:Nt-1
  #   [zeros(Nx, Nu*i+(Nx+Nu)*(i-1)) A B -I
  
  H = sparse([kron(Diagonal(I,Nt-1),[R zeros(Nu,Nx); zeros(Nx,Nu) Q]) zeros((Nx+Nu)*(Nt-1), Nx+Nu); zeros(Nx+Nu,(Nx+Nu)*(Nt-1)) [R zeros(Nu,Nx); zeros(Nx,Nu) Qf]])
  b = zeros(Nt*(Nx+Nu))
  C = sparse([
          [B -I zeros(Nx,(Nt-1)*(Nu+Nx))]; 
          zeros(Nx*(Nt-1),Nu) [kron(Diagonal(I,Nt-1), [A B]) zeros((Nt-1)*Nx,Nx)] + [zeros((Nt-1)*Nx,Nx) kron(Diagonal(I,Nt-1),[zeros(Nx,Nu) Diagonal(-I,Nx)])]
  ])
  Z = kron(Diagonal(I,Nt), [0 0 0 0 0 1 0 0 0]) #Matrix that picks out all x2 (height), for floor constraint
  U = kron(Diagonal(I,Nt), [I zeros(Nu,Nx)]) #Matrix that picks out all u, for thrust constraints
  D = [C; Z; U]
  
  # constraints to keep quadrotor above ground
  lb_floor = zeros(Nt)
  ub_floor = Inf*ones(Nt)

  # thrust constraints
  umin = [-point.max_thrust_x_y, -point.max_thrust_x_y, point.min_thrust_z]-ctrl.Uref[end]
  umax = [point.max_thrust_x_y, point.max_thrust_x_y, point.max_thrust_z]-ctrl.Uref[end]
  lb_thrust = kron(ones(Nt),umin)
  ub_thrust = kron(ones(Nt),umax)

  lb = [zeros(Nx*Nt); lb_floor; lb_thrust]
  ub = [zeros(Nx*Nt); ub_floor; ub_thrust]
  
  Nd = length(ctrl.lb)
  if Nd == Nt*n
      D = C
      ub = zero(ctrl.ub)
      lb = zero(ctrl.lb)
  end
  lb = lb[1:Nd]
  ub = ub[1:Nd]
  
  ctrl.P .= H
  ctrl.A .= D
  ctrl.ub .= ub
  ctrl.lb .= lb
  
  # Initialize the included solver
  #    If you want to use your QP solver, you should write your own
  #    method for this function
  initialize_solver!(ctrl; kwargs...)
  return nothing
end

function updateQP_constrained_point!(ctrl::MPCController, x, time, Q, R, Qf, A)
  # TODO: Implement this method
  
  # SOLUTION
  t = get_k(ctrl, time)
  
  Nt = ctrl.Nmpc-1             # horizon
  Nx = length(ctrl.Xref[1])   # number of states in x̃
  Nu = length(ctrl.Uref[1])    # number of controls
  
  # Update QP problem
  b = ctrl.q
  lb = ctrl.lb
  ub = ctrl.ub
  xeq = ctrl.Xref[end]
  ueq = ctrl.Uref[end]
  N = length(ctrl.Xref)
  for t_h = 1:(Nt-1)
    if (t+t_h) <= N
      # b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*(xref[t+t_h] - xeq)
      # println(size(b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)]))
      # println(size(-Q*(xref[t+t_h]-xeq)))
      # b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*(xref[t+t_h]-xeq)
      b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*(xeq)
      # b[((t_h-1)*(Nx+Nu)).+(1:Nu)] .= -R*(ueq)
    else
      # b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*(xref[end] - xeq)
      # b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*(xref[end]-xeq)
      b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*(xeq)
      # b[((t_h-1)*(Nx+Nu)).+(1:Nu)] .= -R*(-ueq)
    end
  end
  if (t+Nt) <= N
    # b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*(xref[t+Nt] - xeq)
    # println(size(-Qf))
    # b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*(xref[t+Nt]-xeq)
    b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*(xeq)
    # b[((Nt-1)*(Nx+Nu)).+(1:Nu)] .= -R*(ueq)
  else
    # b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*(xref[end] - xeq)
    # b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*(xref[end]-xeq)
    b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*(xeq)
    # b[((Nt-1)*(Nx+Nu)).+(1:Nu)] .= -R*(ueq)
  end
  
  # Update the initial condition
  lb[1:Nx] .= -A*(x)
  ub[1:Nx] .= -A*(x)

  return nothing
end

function simulate_point(model::Point, x0, ctrl, Q, R, Qf, A; tf=ctrl.times[end], dt=2e-2)
  n = 6
  m = 3
  times = range(0, tf, step=dt)
  N = length(times)
  X = [@SVector zeros(n) for k = 1:N] 
  U = [@SVector zeros(m) for k = 1:N-1]
  X[1] = x0

  tstart = time_ns()
  for k = 1:N-1
    if k>200
      ctrl.Xref[end] = [
        5.0, 5.0, 5.0, 
        0.0, 0.0, 0.0
      ]
      A = ForwardDiff.jacobian(_x->point_dynamics_rk4(point, _x, ctrl.Uref[end], dt_point), ctrl.Xref[end])
      B = ForwardDiff.jacobian(_u->point_dynamics_rk4(point, ctrl.Xref[end], _u, dt_point), ctrl.Uref[end])
      buildQP_constrained_point!(point, ctrl, A, B, Q, R, Qf)
    end
    U[k] = get_control_point(ctrl, X[k], times[k], Q, R, Qf, A)
    # u = clamp(U[k], umin, umax)
    X[k+1] = point_dynamics_rk4(model, X[k], U[k], dt)
  end
  tend = time_ns()
  rate = N / (tend - tstart) * 1e9
  println("Controller ran at $rate Hz")
  return X,U,times
end

function simulate_one_step_point(model::Point, x0, ctrl, Q, R, Qf, A, k; tf=ctrl.times[end], dt=2e-2)
  times = range(0, tf, step=dt)

  tstart = time_ns()
  U = get_control_point(ctrl, x0, times[k], Q, R, Qf, A)
  # u = clamp(U[k], umin, umax)
  X = point_dynamics_rk4(model, x0, U, dt)
  tend = time_ns()
  rate = (tend - tstart) * 1e9
  println("Controller ran at $rate Hz")
  return X,U,times
end