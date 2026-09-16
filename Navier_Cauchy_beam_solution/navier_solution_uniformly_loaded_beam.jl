using ModelingToolkit
using MethodOfLines
using DomainSets
using NonlinearSolve

# 1. Define symbolic variables and parameters
# We use @parameters for independent variables (spatial coordinates)
@parameters x y

# We use @variables for the dependent variables (the displacements)
# The (..) signifies they are functions of the independent variables
@variables u(..) v(..)

# 2. Define material parameters
# We can define them as symbolic parameters or as simple constants
E = 1.0e3  # Young's Modulus
ν = 0.3    # Poisson's Ratio

# Calculate Lamé parameters for plane stress
μ = E / (2 * (1 + ν))
λ = (E * ν) / (1 - ν^2)

# 3. Define differential operators
Dx = Differential(x)
Dy = Differential(y)
Dxx = Differential(x)^2
Dyy = Differential(y)^2
Dxy = Differential(x) * Differential(y)

# 4. Write the governing equations
# Divergence: ∇ ⋅ u
div_u = Dx(u(x, y)) + Dy(v(x, y))

# Laplacian: ∇²u and ∇²v
lap_u = Dxx(u(x, y)) + Dyy(u(x, y))
lap_v = Dxx(v(x, y)) + Dyy(v(x, y))

# Navier equation 1 (x-direction)
eq_u = μ * lap_u + (λ + μ) * Dx(div_u) ~ 0

# Navier equation 2 (y-direction)
eq_v = μ * lap_v + (λ + μ) * Dy(div_u) ~ 0

eqs = [eq_u, eq_v]

# 5. Define the domain
# A 2x1 rectangle
L = 2.0
W = 1.0
domains = [x ∈ Interval(0.0, L),
           y ∈ Interval(0.0, W)]

# 6. Define the boundary conditions (BCs)
# This example models a cantilever beam or plate:
# - Left edge (x=0) is fixed (Dirichlet BC)
# - Right edge (x=L) is pulled in the x-direction (Dirichlet BC)
# - Top and bottom edges (y=0, y=W) are free (Neumann BC)
#   (For simplicity, we'll set zero displacement in y on the right
#    and zero shear stress on the free edges)

# Note: The free-surface Neumann BCs are σ_yy = 0 and σ_xy = 0.
# σ_yy = λ(∂u/∂x + ∂v/∂y) + 2μ(∂v/∂y)
# σ_xy = μ(∂u/∂y + ∂v/∂x)
# These are more complex to implement. For a simpler demo,
# we'll use simple Dirichlet and zero-derivative Neumann BCs.

# uniform_load = -10.0 #units are force/length/length


bcs = [
    # Left wall (x=0) is fixed
    u(0, y) ~ 0.0,
    v(0, y) ~ 0.0,

    # Right wall (x=L) is fixed
    u(L, y) ~ 0.0,
#     v(L, y) ~ 0.0, # (restrained in y)
    
    #displace middle of beam 
    v(L, y) ~ -1.0,

    # Top (y=W) and Bottom (y=0) boundaries are free (zero-gradient)
    # This is a simplification of the true "free" stress boundary
        #     λ(∂u/∂x + ∂v/∂y) + 2μ(∂v/∂y)
    Dy(u(x, 0)) ~ 0.0,
    Dy(v(x, 0)) ~ 0.0,
#     λ*(Dx(u(x,W)) + Dy(v(x,W))) + 2 * μ * Dy(v(x,W)) ~ uniform_load
#     Dy(u(x, W)) ~ 0.0,
#     Dy(v(x, W)) ~ 0.0
]

# 7. Create the PDESystem
# Pass the independent and dependent variables
@named pdesys = PDESystem(eqs, bcs, domains, [x, y], [u(x, y), v(x, y)])

# 8. Discretize the system using Method of Lines
# We need to define the grid spacing
dx = 0.1
dy = 0.1

# MOLFiniteDifference handles the discretization
# We pass `nothing` for the time variable to signify a steady-state problem
discretization = MOLFiniteDifference([x => dx, y => dy], nothing)

# discretize converts the PDESystem into a problem type
# For a steady-state problem, this will be a NonlinearProblem
prob = discretize(pdesys, discretization)

# 9. Solve the NonlinearProblem
# We use a solver from NonlinearSolve.jl
println("System discretized. Solving nonlinear problem...")
# sol = NonlinearSolve.solve(prob, NewtonRaphson(), tol = 1e-8)
sol = NonlinearSolve.solve(prob, NewtonRaphson())
println("Solution complete.")

# 10. Extract and plot the solution (Optional)
# The solution `sol` is a flat vector. MethodOfLines provides
# utilities to reconstruct the 2D solution fields.

# Get the discrete grid points
grid_x = sol[x]
grid_y = sol[y]

# Get the solution for u and v
# This returns a 2D array matching the grid
sol_u = sol[u(x, y)]
sol_v = sol[v(x, y)]

println("Solution for u at (x=L, y=W/2): ", sol_u[end, length(grid_y) ÷ 2])

# Optional: Plotting with Plots.jl
# Pkg.add("Plots")
using Plots

heatmap(grid_x, grid_y, sol_u, title="Displacement u(x,y)",
        xlabel="x", ylabel="y", aspect_ratio=:equal)

heatmap(grid_x, grid_y, sol_v, title="Displacement v(x,y)",
        xlabel="x", ylabel="y", aspect_ratio=:equal)