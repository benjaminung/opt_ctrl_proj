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

function get_control(ctrl::MPCController{OSQP.Model}, x, time)
  # Update the QP
  updateQP_constrained!(ctrl, x, time)
  OSQP.update!(ctrl.solver, q=ctrl.q, l=ctrl.lb, u=ctrl.ub)

  # Solve QP
  results = OSQP.solve!(ctrl.solver)
  Δu = results.x[1:4]
  
  k = get_k(ctrl, time)
  return ctrl.Uref[k] + Δu 
end

function buildQP_constrained!(ctrl::MPCController, A,B,Q,R,Qf; kwargs...)
  # TODO: Implement this method to build the QP matrices
  
  # SOLUTION:
  Nt = ctrl.Nmpc-1
  Nx = length(ctrl.Xref[1])    # number of states
  Nu = length(ctrl.Uref[1])    # number of controls

  max_roll = ctrl.max_roll
  min_roll = ctrl.min_roll
  max_pitch = ctrl.max_pitch
  min_pitch = ctrl.min_pitch
  xeq = Xref[end]
  aeq = xeq[7:8]  # actuators
  
  H = sparse([kron(Diagonal(I,Nt-1),[R zeros(Nu,Nx); zeros(Nx,Nu) Q]) zeros((Nx+Nu)*(Nt-1), Nx+Nu); zeros(Nx+Nu,(Nx+Nu)*(Nt-1)) [R zeros(Nu,Nx); zeros(Nx,Nu) Qf]])
  b = zeros(Nt*(Nx+Nu))
  C = sparse([
          [B -I zeros(Nx,(Nt-1)*(Nu+Nx))]; 
          zeros(Nx*(Nt-1),Nu) [kron(Diagonal(I,Nt-1), [A B]) zeros((Nt-1)*Nx,Nx)] + [zeros((Nt-1)*Nx,Nx) kron(Diagonal(I,Nt-1),[zeros(Nx,Nu) Diagonal(-I,Nx)])]
  ])
  Z = kron(Diagonal(I,Nt), [0 0 1 0 0 0 0 0 0 0 0 0]) #Matrix that picks out all x2 (height), for don't go through floor constraint
  Θ = kron(Diagonal(I,Nt), [0 0 0 1 1 1 0 0 0 0 0 0]) #Matrix that picks out the rodrigues param, for pitch and roll constraint
  U = kron(Diagonal(I,Nt), [zeros(Nu,Nx) I]) #Matrix that picks out all u
  D = [C; Z; Θ; U]
  
  umin = ctrl.min_thrust
  umax = ctrl.max_thrust
  lb = [zeros(Nx*Nt); zeros(Nt); -max_roll*ones(Nt); kron(ones(Nt),umin-aeq)]
  ub = [zeros(Nx*Nt); Inf*ones(Nt); max_roll*ones(Nt); kron(ones(Nt),umax-aeq)]
  
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

function updateQP_constrained!(ctrl::MPCController, x, time, xeq)
  # TODO: Implement this method
  
  # SOLUTION
  t = get_k(ctrl, time)
  
  Nt = ctrl.Nmpc-1             # horizon
  Nx = length(ctrl.Xref[1])    # number of states
  Nu = length(ctrl.Uref[1])    # number of controls
  
  # Update QP problem
  b = ctrl.q
  lb = ctrl.lb
  ub = ctrl.ub
  xref = ctrl.Xref
  xeq = Xref[end]
  N = length(ctrl.Xref)
  for t_h = 1:(Nt-1)
    if (t+t_h) <= N
      # b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*(xref[t+t_h] - xeq)
      b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*get_Δx̃(xref[t+t_h],xeq)
    else
      # b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*(xref[end] - xeq)
      b[(Nu+(t_h-1)*(Nx+Nu)).+(1:Nx)] .= -Q*get_Δx̃(xref[end],xeq)
    end
  end
  if (t+Nt) <= N
    # b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*(xref[t+Nt] - xeq)
    b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*get_Δx̃(xref[t+Nt],xeq)
  else
    # b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*(xref[end] - xeq)
    b[(Nu+(Nt-1)*(Nx+Nu)).+(1:Nx)] .= -Qf*get_Δx̃(xref[end],xeq)
  end
  
  # Update the initial condition
  lb[1:Nx] .= -A*get_Δx̃(x,xeq)
  ub[1:Nx] .= -A*get_Δx̃(x,xeq)

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

function simulate(model::Quadrotor, x0, ctrl; tf=ctrl.times[end], dt=1e-2)
    n = 13
    m = 4
    times = range(0, tf, step=dt)
    N = length(times)
    X = [@SVector zeros(n) for k = 1:N] 
    U = [@SVector zeros(m) for k = 1:N-1]
    X[1] = x0

    tstart = time_ns()
    for k = 1:N-1
      U[k] = get_control(ctrl, X[k], times[k])
      # u = clamp(U[k], umin, umax)
      X[k+1] = quad_dynamics_rk4(model, X[k], U[k], dt)
    end
    tend = time_ns()
    rate = N / (tend - tstart) * 1e9
    println("Controller ran at $rate Hz")
    return X,U,times
end