using ForwardDiff
using Plots: scatter3d!
using Makie: scatter!, mesh!, Figure, Axis3, axislegend, contour3d!, contour!
using GLMakie: poly!, contour3d!, contour!, Colorbar, save
using GeometryBasics: Point3f, Mesh

# Set working directory 
cd("C:/Users/chilu/OneDrive - Johns Hopkins/Microsoft VS Code")

H = 12.0            # Height of Bearing (in)
W = 27.0            # Width of Bearing (in)
μ = 0.5             # Poisson's ratio for natural rubber 
a = -2.0            # Displacement (in)
b = -μ*a            # Maximum bulging due to axial displacement (in)

# Generate a mesh grid for a cube (Increase length for finer mesh)
x = range(0, W, length=100)              # u-direction 
y = range(0, W, length=100)              # w-direction
z = range(0, H, length=100)              # v-direction

# Creates a 3D array of grid points representing coordinates of the cube
gP = [Point3f(xi, yi, zi) for xi in x, yi in y, zi in z]

# Defining a function to compute deformations at each grid point
deform(point) = (
    b*sin(π*point[3]/H)*(1/(W/2))*(point[1]-W/2),            # Deformation in the u-direction (x-axis)
    b*sin(π*point[3]/H)*(1/(W/2))*(point[2]-W/2),            # Deformation in the w-direction (y-axis)
    a*point[3]/H                                             # Deformation in the v-direction (z-axis)                     
)

# Apply deformation to each grid point
deformed_points = [Point3f(
    p[1] + deform(p)[1], 
    p[2] + deform(p)[2], 
    p[3] + deform(p)[3]
    ) for p in gP
]

# Visualize the original and deformed grid points
bearing = Figure()
ax = Axis3(bearing[1, 1])
mesh!(ax, [p[1] for p in vec(gP)], [p[2] for p in vec(gP)], [p[3] for p in vec(gP)], color = :blue, alpha = 0.5, label="Original")
mesh!(ax, [p[1] for p in vec(deformed_points)], [p[2] for p in vec(deformed_points)], [p[3] for p in vec(deformed_points)], color = :red, alpha = 0.5, label="Deformed")
axislegend(ax)
display(bearing)
#save("deformed_vs_original_bearing.png", fig)

# Calculate the gradient of the deformation at any point in the domain
strain_u = [Point3f(
      ForwardDiff.gradient(p -> deform(p)[1], point)     # ∇u (x-axis)
)
      for point in gP
]

strain_w = [Point3f(
      ForwardDiff.gradient(p -> deform(p)[2], point)     # ∇w (y-axis)
)
      for point in gP
]

strain_v = [Point3f(
      ForwardDiff.gradient(p -> deform(p)[3], point)     # ∇v (z-axis)
)
      for point in gP
]

# Displacement Gradient Matrix 
disp_grad = [hcat(strain_u[i,j,k], strain_w[i,j,k], strain_v[i,j,k]) for i in 1:length(x), j in 1:length(y), k in 1:length(z)]

# Strain Tensor Matrix (100, 100, 100) 3D array of tensors 
ε_ij = [0.5 * (disp_grad[i,j,k] + disp_grad[i,j,k]') for i in 1:length(x), j in 1:length(y), k in 1:length(z)]

# Extracting the deformed coordinates for plotting
x_def = vec([deformed_points[i,j,k][1] for i in 1:length(x), j in 1:length(y), k in 1:length(z)])
y_def = vec([deformed_points[i,j,k][2] for i in 1:length(x), j in 1:length(y), k in 1:length(z)])
z_def = vec([deformed_points[i,j,k][3] for i in 1:length(x), j in 1:length(y), k in 1:length(z)])

# Extracting the Normal strain components for plotting the strain contours 
ε_xx = vec([ε_ij[i,j,k][1,1] for i in 1:length(x), j in 1:length(y), k in 1:length(z)])               # Normal Strain in u-direction (ε_xx)
ε_yy = vec([ε_ij[i,j,k][2,2] for i in 1:length(x), j in 1:length(y), k in 1:length(z)])               # Normal Strain in w-direction (ε_yy)
ε_zz = vec([ε_ij[i,j,k][3,3] for i in 1:length(x), j in 1:length(y), k in 1:length(z)])               # Normal Strain in v-direction (ε_zz)           

# Extracting the Shear strain components for plotting the strain contours (Assumes ij = ji, so only 3 uniuque shear strains are needed)
ε_xy = vec([ε_ij[i,j,k][1,2] for i in 1:length(x), j in 1:length(y), k in 1:length(z)])               # Shear Strain in u-w plane (ε_xy)
ε_xz = vec([ε_ij[i,j,k][1,3] for i in 1:length(x), j in 1:length(y), k in 1:length(z)])               # Shear Strain in u-v plane (ε_xz)
ε_yz = vec([ε_ij[i,j,k][2,3] for i in 1:length(x), j in 1:length(y), k in 1:length(z)])               # Shear Strain in w-v plane (ε_yz)

# Plots the strain contours for ε_xx
strain_xx = Figure()
ax_xx = Axis3(strain_xx[1,1])
s_xx = scatter!(ax_xx, x_def, y_def, z_def, color=ε_xx, colormap=:viridis, markersize=10, alpha = 0.5)
Colorbar(strain_xx[1,2], s_xx, label="ε_xx")
display(strain_xx)
#save("strain_xx_contours.png", strain_xx)

# Plots the strain contours for ε_yy
strain_yy = Figure()
ax_yy = Axis3(strain_yy[1,1])
s_yy = scatter!(ax_yy, x_def, y_def, z_def, color=ε_yy, colormap=:viridis, markersize=10, alpha = 0.5)
Colorbar(strain_yy[1,2], s_yy, label="ε_yy")
display(strain_yy)
#save("strain_yy_contours.png", strain_yy)

# Plots the strain contours for ε_zz
strain_zz = Figure()
ax_zz = Axis3(strain_zz[1,1])
s_zz = scatter!(ax_zz, x_def, y_def, z_def, color=ε_zz, colormap=:viridis, markersize=10, alpha = 0.5)
Colorbar(strain_zz[1,2], s_zz, label="ε_zz")
display(strain_zz)
#save("strain_zz_contours.png", strain_zz)

# Plots the strain contours for ε_xy
strain_xy = Figure()
ax_xy = Axis3(strain_xy[1,1])
s_xy = scatter!(ax_xy, x_def, y_def, z_def, color=ε_xy, colormap=:viridis, markersize=10, alpha = 0.5)
Colorbar(strain_xy[1,2], s_xy, label="ε_xy")
display(strain_xy)
#save("strain_xy_contours.png", strain_xy)

# Plots the strain contours for ε_yz
strain_yz = Figure()
ax_yz = Axis3(strain_yz[1,1])
s_yz = scatter!(ax_yz, x_def, y_def, z_def, color=ε_yz, colormap=:viridis, markersize=10, alpha = 0.5)
Colorbar(strain_yz[1,2], s_yz, label="ε_yz")
display(strain_yz)
#save("strain_yz_contours.png", strain_yz)

# Plots the strain contours for ε_xz
strain_xz = Figure()
ax_xz = Axis3(strain_xz[1,1])
s_xz = scatter!(ax_xz, x_def, y_def, z_def, color=ε_xz, colormap=:viridis, markersize=10, alpha = 0.5)
Colorbar(strain_xz[1,2], s_xz, label="ε_xz")
display(strain_xz)
#save("strain_xz_contours.png", strain_xz)

#=
# Manually computes the displacement gradient and strain tensor at a specific point 
point = [W, W, H/2]                       # [X, Y, Z] or [u, w, v]

# Displacement Gradient Matrix
disp_grad = zeros(3,3)
for i in 1:3
      for j in 1:3
            disp_grad[i,j] = [∇u[i], ∇w[i], ∇v[i]][j]
      end
end
disp_grad = disp_grad'

# Strain Tensor Matrix 
ε_ij = zeros(3,3)
for i in 1:3
      for j in 1:3
            ε_ij[i,j] = 0.5*(disp_grad[i,j] + disp_grad[j,i])
      end
end
=#

#=
# Check the vectors and matrices for correctness
print("∇u is $∇u\n")
print("∇w is $∇w\n")
print("∇v is $∇v\n")
display(disp_grad)
display(ε_ij)
=#


