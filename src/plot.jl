using Plots

Y_red = [X_red_sim_const_vel[i][2] for i=1:857]
Z_red = [X_red_sim_const_vel[i][3] for i=1:857]
Y_blue = [X_blue_const_vel[i][2] for i=1:857]
Z_blue = [X_blue_const_vel[i][3] for i=1:857]
graph_const_velocity = plot(Y_blue, Z_blue, c=:blue, xlabel="Y positon [m]", ylabel="Z position [m]", xlims=(-5, 30), ylims=(10,25), title="Red Const Vel", labels="Blue Traj")
plot!(Y_red, Z_red, c=:red, labels="Red Traj")
savefig(graph_const_velocity, "plot_const_vel.png")

Y_red = [X_red_sim_zigz[i][2] for i=1:857]
Z_red = [X_red_sim_zigz[i][3] for i=1:857]
Y_blue = [X_blue_zigz[i][2] for i=1:857]
Z_blue = [X_blue_zigz[i][3] for i=1:857]
graph_zigz = plot(Y_blue, Z_blue, c=:blue, xlabel="Y positon [m]", ylabel="Z position [m]", xlims=(-5, 30), ylims=(10,25), title="Red Zig Z", labels="Blue Traj")
plot!(Y_red, Z_red, c=:red, labels="Red Traj")
savefig(graph_zigz, "plot_zigz.png")

Y_red = [X_red_sim_zigy[i][2] for i=1:837]
Z_red = [X_red_sim_zigy[i][3] for i=1:837]
Y_blue = [X_blue_zigy[i][2] for i=1:837]
Z_blue = [X_blue_zigy[i][3] for i=1:837]
graph_zigy = plot(Y_blue, Z_blue, c=:blue, xlabel="Y positon [m]", ylabel="Z position [m]", xlims=(-5, 30), ylims=(10,25), title="Red Zig Y", labels="Blue Traj")
plot!(Y_red, Z_red, c=:red, labels="Red Traj")
savefig(graph_zigy, "plot_zigy.png")

Y_red = [X_red_sim_zigx[i][2] for i=1:866]
X_red = [X_red_sim_zigx[i][1] for i=1:866]
Y_blue = [X_blue_zigx[i][2] for i=1:866]
X_blue = [X_blue_zigx[i][1] for i=1:866]
graph_zigx = plot(X_blue, Y_blue, c=:blue, xlabel="X positon [m]", ylabel="Y position [m]", xlims=(-5, 5), ylims=(10,25), title="Red Zig X", labels="Blue Traj")
plot!(X_red, Y_red, c=:red, labels="Red Traj")
savefig(graph_zigx, "plot_zigx.png")