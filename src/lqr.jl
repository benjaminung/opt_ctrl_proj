struct LQRController
  K::Vector{Matrix{Float64}}   # feedback gains ((m,n),N-1)
  times::Vector{Float64}       # times          (N,)
  Xref::Vector{Vector{Float64}}         # equilibrium states
  Uref::Vector{Vector{Float64}}         # equilibrium controls
end

function riccati(A,B,Q,R,Qf)
  n,m = size(B[1])
  N = length(A)+1
  P = [zeros(n,n) for k = 1:N]
  K = [zeros(m,n) for k = 1:N-1]
  
  P[end] .= Qf
  for k = reverse(1:N-1) 
      K[k] .= (R + B[k]'P[k+1]*B[k])\(B[k]'P[k+1]*A[k])
      P[k] .= Q + A[k]'P[k+1]*A[k] - A[k]'P[k+1]*B[k]*K[k]
  end
  
  return K,P
end

function get_Ã_B̃(quad, xeq0, xeq1, ueq, Q̃, R, dt)
  E_xeq0 = E(xeq0[4:7])
  E_xeq1 = E(xeq1[4:7])
  A = ForwardDiff.jacobian(_x->quad_dynamics_rk4(quad, _x, ueq, dt), xeq0)
  Ã = E_xeq1' * A * E_xeq0
  B = ForwardDiff.jacobian(_u->quad_dynamics_rk4(quad, xeq0, _u, dt), ueq)
  B̃ = E_xeq1' * B
  return Ã, B̃
end

function controller(x, xeq, ueq, K) 
  qeq = xeq[4:7]
  q = x[4:7]
  ϕ = qtorp(Lmult(qeq)'*q)
  
  Δx̃ = [x[1:3]-xeq[1:3]; ϕ; x[8:10]-xeq[8:10]; x[11:13]-xeq[11:13]]
  
  u = ueq - K*Δx̃
end

function get_Δx̃(x, xeq)
  qeq = xeq[4:7]
  q = x[4:7]
  ϕ = qtorp(Lmult(qeq)'*q)
  
  Δx̃ = [x[1:3]-xeq[1:3]; ϕ; x[8:10]-xeq[8:10]; x[11:13]-xeq[11:13]]
  return Δx̃
end

get_k(controller::LQRController, t) = searchsortedlast(controller.times, t)

function get_control(controller::LQRController, x, k, min_u=-Inf, max_u=Inf)
  u = controller.Uref[k] - controller.K[k]*get_Δx̃(x, controller.Xref[k])
  u = min.(u, max_u)
  u = max.(u, min_u)
  return u
end