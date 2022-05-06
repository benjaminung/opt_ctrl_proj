using StaticArrays
using ControlSystems
using OSQP

struct LQRController
  K::Vector{Matrix{Float64}}
  Xref::Vector{Vector{Float64}}
  Uref::Vector{Vector{Float64}}
  times::Vector{Float64}
end

function riccati(A,B,Q,R,Qf,N)
  # initialize the output
  n,m = size(B)
  P = [zeros(n,n) for k = 1:N]
  K = [zeros(m,n) for k = 1:N-1]

  # TODO: implement the Riccati recursion
  P[end] .= Qf
  for k = reverse(1:N-1) 
      K[k] .= (R + B'P[k+1]*B)\(B'P[k+1]*A)
      P[k] .= Q + A'P[k+1]*A - A'P[k+1]*B*K[k]
  end
  
  # return the feedback gains and ctg matrices
  return K,P
end

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