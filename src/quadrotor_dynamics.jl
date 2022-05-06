using LinearAlgebra

struct Quadrotor
  m::Float64  # mass [kg]
  ℓ::Float64  # distance between motors [meters]
  J::SMatrix{3,3,Float64}  # inerita [kg m^2]
  g::Float64  # gravity acceleration [m/s^2]
  kt::Float64  # force constant
  km::Float64  # torque constant
  max_thrust::Float64
  min_thrust::Float64
end

function Quadrotor()
  m = 0.5
  ℓ = 0.1750
  J = Diagonal([0.0023, 0.0023, 0.004])
  g = 9.81
  kt = 1.0
  km = 0.0245
  max_thrust = 0.6*m*g
  min_thrust = 0.1*m*g
  Quadrotor(m,ℓ,J,g,kt,km,max_thrust,min_thrust)
end

function quad_dynamics(quad::Quadrotor,x,u)
  H = [zeros(1,3); I]
  m = quad.m
  J = quad.J
  ℓ = quad.ℓ
  kt = quad.kt
  km = quad.km
  g = quad.g

  r = x[1:3]
  q = x[4:7]
  v = x[8:10]
  ω = x[11:13]
  Q = qtoQ(q)
  
  ṙ = Q*v
  q̇ = 0.5*Lmult(q)*H*ω
  
  # v̇ = Q'*[0; 0; -g] + (1/m)*[zeros(2,4); kt*ones(1,4)]*u - hat(ω)*v # body velocity
  v̇ = [0; 0; -g] + Q*(1/m)*[zeros(2,4); kt*ones(1,4)]*u - Q*hat(ω)*v # world velocity
  
  ω̇ = J\(-hat(ω)*J*ω + [0 ℓ*kt 0 -ℓ*kt; -ℓ*kt 0 ℓ*kt 0; km -km km -km]*u)
  
  return [ṙ; q̇; v̇; ω̇]
end

function quad_dynamics_rk4(quad,x,u,h)
  #RK4 integration with zero-order hold on u
  f1 = quad_dynamics(quad, x, u)
  f2 = quad_dynamics(quad, x + 0.5*h*f1, u)
  f3 = quad_dynamics(quad, x + 0.5*h*f2, u)
  f4 = quad_dynamics(quad, x + h*f3, u)
  xn = x + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)
  xn[4:7] .= xn[4:7]/norm(xn[4:7]) #re-normalize quaternion
  return xn
end