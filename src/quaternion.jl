using LinearAlgebra
using BlockDiagonals

#Quaternion stuff
function hat(v)
  return [0 -v[3] v[2];
          v[3] 0 -v[1];
          -v[2] v[1] 0]
end

function Lmult(q)
  s = q[1]
  v = q[2:4]
  L = [s    -v';
       v  s*I+hat(v)]
  return L
end

function qtoQ(q)
  T = Diagonal([1; -ones(3)])
  H = [zeros(1,3); I]
  return H'*T*Lmult(q)*T*Lmult(q)*H
end

function G(q)
  H = [zeros(1,3); I]
  G = Lmult(q)*H
end

function rptoq(ϕ)
  (1/sqrt(1+ϕ'*ϕ))*[1; ϕ]
end

function x_to_x̃(x)
  [x[1:3]; qtorp(x[4:7]); x[8:13]]
end

function body_to_world_vel(x)
  Q = qtoQ(x[4:7])
  v = x[8:10]
  return Q*v
end

function qtorp(q)
  q[2:4]/q[1]
end

function E(q)
  E = BlockDiagonal([1.0*I(3), G(q), 1.0*I(6)])
end