using MeshCat
using CoordinateTransformations
using Rotations 
using Colors
using GeometryBasics: HyperRectangle, Vec
using StaticArrays

function set_mesh!(vis)
    urdf_folder = @__DIR__
    obj = joinpath(urdf_folder, "../quadrotor_scaled.obj")
    robot_obj = MeshFileGeometry(obj)
    defend_obj = HyperRectangle(Vec(-0.25, -0.25, 0.0), Vec(0.5, 0.5, 0.5))

    red_mat = MeshPhongMaterial(color=colorant"red")
    blue_mat = MeshPhongMaterial(color=colorant"blue")
    green_mat = MeshPhongMaterial(color=colorant"green")
    setobject!(vis["red"]["geom"], robot_obj, red_mat)
    setobject!(vis["blue"]["geom"], robot_obj, blue_mat)
    setobject!(vis["defend"]["geom"], defend_obj, green_mat)
end

function visualize!(vis, tf::Real, X_blue, X_red)
    fps = Int(round((length(X_blue)-1)/tf))
    # fps = 50
    anim = MeshCat.Animation(fps)
    for (k,x) in enumerate(X_blue)
        # if k%2 == 0
        #     continue
        # end
        atframe(anim, k) do
            x_blue = X_blue[k]
            x_red = X_red[k]
            visualize!(vis, "blue", SVector{13}(x_blue))
            visualize!(vis, "red", SVector{13}(x_red))
        end
    end
    setanimation!(vis, anim)
end

function visualize!(vis, vis_obj_name::String, x::StaticArray)
    settransform!(vis[vis_obj_name], compose(Translation(x[1]/10, x[2]/10, x[3]/10), LinearMap(UnitQuaternion(x[4:7]))))
    # settransform!(vis[vis_obj_name], Translation(x[1]/10, x[2]/10, x[3]/10))
end

function initialize_visualizer()
    vis = Visualizer()
    set_mesh!(vis)
    return vis
end