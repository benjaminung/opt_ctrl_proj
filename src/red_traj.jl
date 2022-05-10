@enum RedTraj begin
  const_vel
  zig_zag_z
  zig_zag_x
  zig_zag_y
end

function red_traj(red::Quadrotor, times, redTraj::RedTraj=const_vel, frequency=0.2, amplitude=5.0)
  dt = times[2]-times[1]

  red_cost_Q̃ = Diagonal([10, 10, 10, 1, 1, 1, 10, 10, 10, 1, 1, 1])
  red_cost_R = Diagonal([0.1, 0.1, 0.1, 0.1])

  uhover = [0.0, 0.0, 0.0, 0.0] .+ red.m*red.g/4.0
  u0_red = copy(uhover)
  Uref_red = [copy(u0_red) for i=1:length(times)-1]
  
  pos0 = [0.0, 100.0, 100.0]
  local Xref_red
  if redTraj == const_vel
    Xref_red = red_Xref_constant_vel(times, pos0)
  elseif redTraj == zig_zag_z
    Xref_red = red_Xref_zig_zag_z(times, pos0, frequency, amplitude)
  elseif redTraj == zig_zag_x
    Xref_red = red_Xref_zig_zag_x(times, pos0, frequency, amplitude)
  elseif redTraj == zig_zag_y
    Xref_red = red_Xref_zig_zag_y(times, pos0, frequency, amplitude)
  end

  A = [zeros(12,12) for i=1:length(times)-1]
  B = [zeros(12,4) for i=1:length(times)-1]

  for i=1:length(times)-1
    A[i], B[i] = get_Ã_B̃(red, Xref_red[i], Xref_red[i+1], Uref_red[i], red_cost_Q̃, red_cost_R, dt)
  end
  red_K, = riccati(A, B, red_cost_Q̃, red_cost_R, red_cost_Q̃)

  red_lqr = LQRController(red_K, times, Xref_red, Uref_red)

  x0_red = Xref_red[1]
  X_red = [copy(x0_red) for i=1:length(times)]
  U_red = [copy(uhover) for i=1:length(times)-1]

  for i=1:length(times)-1
    U_red[i] = get_control(red_lqr, X_red[i], i)
    X_red[i+1] = quad_dynamics_rk4(red, X_red[i], U_red[i], times[i+1]-times[i])
  end

  return X_red
end

function red_Xref_constant_vel(times, initial_position)
  dt = times[2]-times[1]

  quat0 = [1.0,  0.0,  0.0,  0.0]
  velocity0 = -initial_position/times[end]
  ω0 = [0.0,  0.0,  0.0]

  x0_red = [initial_position; quat0; velocity0; ω0]  
  Xref_red = [copy(x0_red) for i=1:length(times)]
  for i=1:length(times)-1
    xref = copy(Xref_red[i])
    Xref_red[i+1][1:3] = xref[1:3] + Xref_red[i][8:10] * dt
  end
  
  return Xref_red
end

function red_Xref_zig_zag_z(times, initial_position, frequency=0.2, amplitude=5.0)
  dt = times[2]-times[1]

  quat0 = [1.0,  0.0,  0.0,  0.0]
  velocity0 = -initial_position/times[end]
  ω0 = [0.0,  0.0,  0.0]

  x0_red = [initial_position; quat0; velocity0; ω0]
  Xref_red = [copy(x0_red) for i=1:length(times)]
  for i=1:length(times)-1
    Xref_red[i][10] = copy(Xref_red[i][10]) + amplitude*cos(2*pi*times[i]*frequency + pi)
    xref = copy(Xref_red[i])
    Xref_red[i+1][1:3] = xref[1:3] + Xref_red[i][8:10] * dt
  end

  return Xref_red
end

function red_Xref_zig_zag_y(times, initial_position, frequency=0.2, amplitude=5.0)
  dt = times[2]-times[1]

  quat0 = [1.0,  0.0,  0.0,  0.0]
  velocity0 = -initial_position/times[end]
  ω0 = [0.0,  0.0,  0.0]

  x0_red = [initial_position; quat0; velocity0; ω0]
  Xref_red = [copy(x0_red) for i=1:length(times)]
  for i=1:length(times)-1
    Xref_red[i][9] = copy(Xref_red[i][9]) + amplitude*cos(2*pi*times[i]*frequency + pi)
    xref = copy(Xref_red[i])
    Xref_red[i+1][1:3] = xref[1:3] + Xref_red[i][8:10] * dt
  end

  return Xref_red
end

function red_Xref_zig_zag_x(times, initial_position, frequency=0.2, amplitude=5.0)
  dt = times[2]-times[1]

  quat0 = [1.0,  0.0,  0.0,  0.0]
  velocity0 = -initial_position/times[end]
  ω0 = [0.0,  0.0,  0.0]

  x0_red = [initial_position; quat0; velocity0; ω0]
  Xref_red = [copy(x0_red) for i=1:length(times)]
  for i=1:length(times)-1
    Xref_red[i][8] = copy(Xref_red[i][8]) + amplitude*cos(2*pi*times[i]*frequency + pi)
    xref = copy(Xref_red[i])
    Xref_red[i+1][1:3] = xref[1:3] + Xref_red[i][8:10] * dt
  end

  return Xref_red
end