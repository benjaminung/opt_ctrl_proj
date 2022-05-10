using SparseArrays
using OSQP

struct MPCController{S}
  P::SparseMatrixCSC{Float64,Int}
  q::Vector{Float64}
  A::SparseMatrixCSC{Float64,Int}
  lb::Vector{Float64}
  ub::Vector{Float64}
  Nmpc::Int
  solver::S
  Xref::Vector{Vector{Float64}}
  Uref::Vector{Vector{Float64}}
  times::Vector{Float64}
end

function OSQPController(n::Integer, m::Integer, N::Integer, Nref::Integer=N, Nd::Integer=(N-1)*n)
  Np = (N-1)*(n+m)   # number of primals
  P = spzeros(Np,Np)
  q = zeros(Np)
  A = spzeros(Nd,Np)
  lb = zeros(Nd)
  ub = zeros(Nd)
  Xref = [zeros(n) for k = 1:Nref]
  Uref = [zeros(m) for k = 1:Nref]
  tref = zeros(Nref)
  solver = OSQP.Model()
  MPCController{OSQP.Model}(P,q, A,lb,ub, N, solver, Xref, Uref, tref)
end

function initialize_solver!(ctrl::MPCController{OSQP.Model}; tol=1e-6, verbose=false)
  OSQP.setup!(ctrl.solver, P=ctrl.P, q=ctrl.q, A=ctrl.A, l=ctrl.lb, u=ctrl.ub, 
      verbose=verbose, eps_rel=tol, eps_abs=tol, polish=1)
end

function get_control(ctrl::MPCController{OSQP.Model}, x, time, Q, Qf, A)
  # Update the QP
  updateQP_constrained!(ctrl, x, time, Q, Qf, A)
  OSQP.update!(ctrl.solver, q=ctrl.q, l=ctrl.lb, u=ctrl.ub)

  # Solve QP
  results = OSQP.solve!(ctrl.solver)
  Δu = results.x[1:4]
  # println(Δu)
  
  k = get_k(ctrl, time)
  u = ctrl.Uref[end] + Δu 
  if maximum(isnan.(u))
    u = ctrl.Uref[end]
  end
  return u
end

function buildQP_constrained!(quad::Quadrotor, ctrl::MPCController, A,B,Q,R,Qf; kwargs...)
  Nt = ctrl.Nmpc-1
  Nx = length(ctrl.Xref[1])-1    # number of states in x̃
  Nu = length(ctrl.Uref[1])    # number of controls
  
  H = sparse([kron(Diagonal(I,Nt-1),[R zeros(Nu,Nx); zeros(Nx,Nu) Q]) zeros((Nx+Nu)*(Nt-1), Nx+Nu); zeros(Nx+Nu,(Nx+Nu)*(Nt-1)) [R zeros(Nu,Nx); zeros(Nx,Nu) Qf]])
  b = zeros(Nt*(Nx+Nu))
  C = sparse([
          [B -I zeros(Nx,(Nt-1)*(Nu+Nx))]; 
          zeros(Nx*(Nt-1),Nu) [kron(Diagonal(I,Nt-1), [A B]) zeros((Nt-1)*Nx,Nx)] + [zeros((Nt-1)*Nx,Nx) kron(Diagonal(I,Nt-1),[zeros(Nx,Nu) Diagonal(-I,Nx)])]
  ])
  Z = kron(Diagonal(I,Nt), [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]) #Matrix that picks out all x2 (height), for floor constraint
  Θ = kron(Diagonal(I,Nt), [zeros(3, Nu + 3) I zeros(3, 6)]) #Matrix that picks out the rodrigues param, for constraint to keep all euler rotation angles under 90 degrees
  U = kron(Diagonal(I,Nt), [I zeros(Nu,Nx)]) #Matrix that picks out all u, for thrust constraints
  D = [C; Z; Θ; U]
  
  # constraints to keep quadrotor above ground
  lb_floor = zeros(Nt)
  ub_floor = Inf*ones(Nt)

  # each part of rp needs to be between -1 and 1 to keep quadrotor euler angles between -90 and 90 degrees
  lb_angle = kron(ones(Nt), -ones(3))
  ub_angle = kron(ones(Nt), ones(3))

  # thrust constraints
  umin = [copy(quad.min_thrust) for i=1:Nu] - ctrl.Uref[end]
  umax = [copy(quad.max_thrust) for i=1:Nu] - ctrl.Uref[end]
  lb_thrust = kron(ones(Nt),umin)
  ub_thrust = kron(ones(Nt),umax)

  lb = [zeros(Nx*Nt); lb_floor; lb_angle; lb_thrust]
  ub = [zeros(Nx*Nt); ub_floor; ub_angle; ub_thrust]
  
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
  
  # Initialize the solver
  initialize_solver!(ctrl; kwargs...)
  return nothing
end

function updateQP_constrained!(ctrl::MPCController, x, time, Q, Qf, A)
  t = get_k(ctrl, time)
  
  Nt = ctrl.Nmpc-1             # horizon
  Nx = length(ctrl.Xref[1]) - 1   # number of states in x̃
  Nu = length(ctrl.Uref[1])    # number of controls
  
  # Update QP problem
  b = ctrl.q
  lb = ctrl.lb
  ub = ctrl.ub
  xref = ctrl.Xref
  xeq = xref[end]
  N = length(ctrl.Xref)
  for t_h = 1:(Nt-1)
    if (t+t_h) <= N
      # b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*(xref[t+t_h] - xeq)
      b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*x_to_x̃(xeq)
    else
      # b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*(xref[end] - xeq)
      b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*x_to_x̃(xeq)
    end
  end
  if (t+Nt) <= N
    # b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*(xref[t+Nt] - xeq)
    b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*x_to_x̃(xeq)
  else
    # b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*(xref[end] - xeq)
    b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*x_to_x̃(xeq)
  end
  
  # Update the initial condition
  lb[1:Nx] .= -A*x_to_x̃(x)
  ub[1:Nx] .= -A*x_to_x̃(x)

  return nothing
end

get_k(controller::MPCController, t) = searchsortedlast(controller.times, t)

function update_xref!(ctrl::MPCController, x, time, dt, xeq, time_eq)
  t = get_k(ctrl, time)
  # println(t)
  teq = get_k(ctrl, time_eq)
  # println(teq)
  vel = (copy(xeq[1:3]) - copy(x[1:3]))/(time_eq-time)
  # println(vel)
  for k=t:teq - 1
    ctrl.Xref[k][8:10] = copy(vel)
    ctrl.Xref[k+1][1:3] = copy(ctrl.Xref[k][1:3]) + copy(ctrl.Xref[k][8:10])*dt
  end
  if teq < length(ctrl.times)
    ctrl.Xref[teq][8:10] = [0.0, 0.0, 0.0]
    for k=teq:length(ctrl.times)-1
      ctrl.Xref[k+1][1:3] = copy(xeq[1:3])
      ctrl.Xref[k+1][8:10] = [0.0, 0.0, 0.0]
    end
  end
end

function simulate(model::Quadrotor, x0, ctrl, Q, R, Qf, A, redTraj; tf=ctrl.times[end], dt=2e-2)
  n = 13
  m = 4
  times = range(0, tf, step=dt)
  N = length(times)
  X = [@SVector zeros(n) for k = 1:N] 
  U = [@SVector zeros(m) for k = 1:N-1]
  X[1] = x0

  tstart = time_ns()
  local intercept_k = 0
  for k = 1:N-1
    int_pos, int_vel = get_intercept_point_vel(redTraj[k], ctrl.Xref[end][4:7])
    ctrl.Xref[end][1:3] = int_pos
    # ctrl.Xref[end][8:10] = int_vel
    A, B = get_Ã_B̃(model, ctrl.Xref[end], ctrl.Xref[end], ctrl.Uref[end], Q, R, dt)
    buildQP_constrained!(model, ctrl, A, B, Q, R, Q)
    U[k] = get_control(ctrl, X[k], times[k], Q, Qf, A)
    X[k+1] = quad_dynamics_rk4(model, X[k], U[k], dt)
    intercepted = interception_check(X[k+1], redTraj[k])
    if intercepted
      intercept_k = k+1
      break 
    end
  end
  tend = time_ns()
  rate = N / (tend - tstart) * 1e9
  println("Controller ran at $rate Hz")
  if intercept_k > 0
    println("INTERCEPTED @ t=", times[intercept_k])
    return X[1:intercept_k], U[1:intercept_k-1], times[1:intercept_k], redTraj[1:intercept_k]
  end
  return X,U,times,redTraj
end

function simulate_one_step(model::Quadrotor, x0, ctrl, Q, Qf, A, k; tf=ctrl.times[end], dt=2e-2)
  times = range(0, tf, step=dt)

  tstart = time_ns()
  U = get_control(ctrl, x0, times[k], Q, Qf, A)
  # u = clamp(U[k], umin, umax)
  X = quad_dynamics_rk4(model, x0, U, dt)
  tend = time_ns()
  rate = (tend - tstart) * 1e9
  println("Controller ran at $rate Hz")
  return X,U,times
end

function get_intercept_point_vel(x_red, quateq_blue, dt=0.02, distance=20.0)
  red_pos = copy(x_red[1:3])
  red_world_vel = body_to_world_vel(x_red[8:10], x_red[4:7])
  if norm(red_pos) > distance
    red_pos_unit = red_pos/norm(red_pos)
    intercept_point = red_pos_unit*distance
    intercept_vel = red_world_vel - red_pos_unit*dot(red_world_vel, red_pos_unit)
    # intercept_vel = intercept_vel * distance/norm(red_pos)
    intercept_point = intercept_point + intercept_vel*dt*50 
    return intercept_point, intercept_vel
  else
    return red_pos, [0, 0, 0]
  end
end

function interception_check(x_blue, x_red, tolerance=2.0)
  dist = norm(x_blue[1:3]-x_red[1:3])
  if dist < tolerance
    return true
  end
  return false
end
  