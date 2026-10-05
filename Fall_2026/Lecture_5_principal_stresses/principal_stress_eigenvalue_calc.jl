# Run from this folder with:  julia --project=. principal_stress_eigenvalue_calc.jl
# or in VS Code, activate this folder's environment first.

using LinearAlgebra
# using WGLMakie

#pure shear example, only σ12 = σ21 = τ is nonzero, units are MPa
τ = 10.0

#tensor
σij = [0.0   τ     0.0
       τ     0.0   0.0
       0.0   0.0   0.0]

#principle stresses, expect -τ, 0, +τ 
σ = eigvals(σij)

#columns in matrix are eigenvectors, expect principal directions at 45° to x1 and x2 
σ_directions = eigvecs(σij)


############################################################################
# plot the stress state and the principal directions on a 3D cube
############################################################################

#cube centered at the origin with side length L 
L = 1.0
cube = Rect3f(Point3f(-L/2, -L/2, -L/2), Vec3f(L, L, L))

fig = Figure(size = (900, 800))
ax = Axis3(fig[1, 1],
           aspect = :data, limits = (-1.3L, 1.3L, -1.3L, 1.3L, -1.3L, 1.3L),
           azimuth = 0.3π, elevation = 0.2π,
           xlabel = "x₁", ylabel = "x₂", zlabel = "x₃",
           title = "pure shear, σ₁₂ = σ₂₁ = τ = $τ MPa")

#translucent cube with black edges
mesh!(ax, cube, color = (:lightsteelblue, 0.25), transparency = true)
wireframe!(ax, cube, color = :black, linewidth = 1.5)

#traction vectors t = σij ⋅ n acting on each of the six cube faces 
#arrows are scaled so the largest stress has length 0.4 L
σ_max = maximum(abs.(σ))
scale = 0.4L / σ_max

normals = [Vec3f(1,0,0), Vec3f(-1,0,0), Vec3f(0,1,0), Vec3f(0,-1,0), Vec3f(0,0,1), Vec3f(0,0,-1)]
face_centers = Point3f[]
tractions = Vec3f[]
for n in normals
    t = σij * n            #Cauchy traction on the face with outward normal n
    if norm(t) > 1e-9      #skip faces with zero traction
        push!(face_centers, Point3f(n * L/2))
        push!(tractions, Vec3f(t * scale))
    end
end
arrows3d!(ax, face_centers, tractions, color = :darkorange, 
          shaftradius = 0.03, tipradius = 0.07, tiplength = 0.12, markerscale = 1.0)

#principal stress directions through the cube center 
#red = tension, blue = compression, gray = zero principal stress 
#eigenvector sign is arbitrary so draw arrows in both ±n directions 
principal_colors = [s > 1e-9 ? :red : s < -1e-9 ? :blue : :gray for s in σ]
for k in 1:3
    n = normalize(σ_directions[:, k])
    for sgn in (1, -1)
        arrows3d!(ax, [Point3f(0, 0, 0)], [Vec3f(sgn * 0.9L * n)], color = principal_colors[k],
                  shaftradius = 0.025, tipradius = 0.06, tiplength = 0.12, markerscale = 1.0)
    end
    text!(ax, Point3f(1.05L * n), text = "σ = $(round(σ[k], digits=2)) MPa", 
          color = principal_colors[k], fontsize = 16, align = (:center, :center))
end

#dotted lines showing the x₁, x₂, x₃ reference axes 
for (e, label) in zip((Vec3f(1,0,0), Vec3f(0,1,0), Vec3f(0,0,1)), ("x₁", "x₂", "x₃"))
    lines!(ax, [Point3f(0,0,0), Point3f(0.8L * e)], color = :black, linestyle = :dot, linewidth = 2)
    text!(ax, Point3f(0.9L * e), text = label, fontsize = 18, align = (:center, :center))
end

display(fig)

#WGLMakie handles mouse rotation of Axis3 on the Julia side, so if this file is run as a 
#script (julia --project=. file.jl) Julia would exit right after display and the browser 
#page would freeze. Keep the process alive until the user presses Enter. 
#Not needed when running in the REPL or VS Code, where Julia stays open.
if !isinteractive()
    println("figure is live in the browser, press Enter here to exit")
    readline()
end
