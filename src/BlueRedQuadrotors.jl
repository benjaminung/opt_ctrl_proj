using RobotZoo


struct BlueRedQuadrotors
  model_blue::RobotZoo.Quadrotor
  model_red::RobotZoo.Quadrotor
  knotpoints_count::Int
  times::Vector{Float64}
  X_blue::Vector{Vector{Float64}}
  U_blue::Vector{Vector{Float64}}
  X_red::Vector{Vector{Float64}}
  U_red::Vector{Vector{Float64}}
end