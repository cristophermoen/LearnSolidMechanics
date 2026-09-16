# Install packages if needed
# using Pkg
# Pkg.add(["ModelingToolkit", "NeuralPDE", "Lux", "Optimization", "OptimizationOptimisers", "DomainSets"])

using ModelingToolkit
using NeuralPDE
using Lux, Optimization, OptimizationOptimisers
using DomainSets
using Plots

# Suppress annoying warnings if any
# import ModelingToolkit: Interval, inf, sup


@parameters x y
@variables u(..) v(..)

# Use `(..)` to signify that u and v are functions of variables defined later
# In this case, they will be functions of x and y.

# Define derivative operators
Dx = Differential(x)
Dy = Differential(y)

# --- Physical Parameters ---
# We will substitute values later, but define them symbolically for now
@parameters E ν q L h

# --- Material Constants ---
# Plane stress material matrix components
C11 = E / (1 - ν^2)
C12 = ν * E / (1 - ν^2)
C66 = E / (2 * (1 + ν)) # This is the shear modulus, G

# --- Kinematics (Strains) ---
ϵ_x = Dx(u(x, y))
ϵ_y = Dy(v(x, y))
γ_xy = Dy(u(x, y)) + Dx(v(x, y))

# --- Constitutive (Stresses) ---
σ_x = C11 * ϵ_x + C12 * ϵ_y
σ_y = C12 * ϵ_x + C11 * ϵ_y
τ_xy = C66 * γ_xy

# --- Equilibrium (PDEs) ---
eq1 = Dx(σ_x) + Dy(τ_xy) ~ 0
eq2 = Dx(τ_xy) + Dy(σ_y) ~ 0

eqs = [eq1, eq2]


# Domain: A rectangle
domains = [x ∈ Interval(0.0, L),
           y ∈ Interval(-h/2, h/2)]


# Boundary Conditions
bcs = [
    # Fixed end at x = 0 (Dirichlet BCs)
    u(0, y) ~ 0,
    v(0, y) ~ 0,

    # Free end at x = L (Neumann BCs, zero stress)
    σ_x(L, y) ~ 0, # Note: We use the symbolic stress σ_x
    τ_xy(L, y) ~ 0,

    # Free bottom at y = -h/2 (Neumann BCs, zero stress)
    σ_y(x, -h/2) ~ 0,
    τ_xy(x, -h/2) ~ 0,

    # Loaded top at y = h/2 (Neumann BCs)
    σ_y(x, h/2) ~ -q, # Uniform downward load
    τ_xy(x, h/2) ~ 0  # No shear load
]


@named pde_system = PDESystem(eqs, bcs, domains, [x, y], [u(x, y), v(x, y)], [E, ν, q, L, h])


# Number of neurons in each hidden layer
inner_layers = 16

# Use Lux.jl to define the network structure
chain = Lux.Chain(Lux.Dense(2, inner_layers, Lux.tanh),  # Input (x,y) -> 16
                  Lux.Dense(inner_layers, inner_layers, Lux.tanh), # 16 -> 16
                  Lux.Dense(inner_layers, 2)) # 16 -> 2 (u, v)


                  # Discretization strategy
discretization = PhysicsInformedNN(chain, QuadratureTraining())

# Create the problem
prob = discretize(pde_system, discretization)


# Assign concrete values
param_values = [
    E => 1000.0,  # Young's Modulus
    ν => 0.3,     # Poisson's Ratio
    q => 1.0,     # Load magnitude
    L => 10.0,    # Beam length
    h => 2.0      # Beam height
]

# Create a mapping for the `discretize` function
# This is a bit of a workaround to map symbolic params to their values
prob = Base.remake(prob, p = param_values)

# Callback function to observe training progress
callback = function (p, l)
    println("Current loss: $l")
    return false
end

# --- Training ---
# Use a strong optimizer
optimizer = OptimizationOptimisers.Adam(0.01)
# Train for 500 iterations
res = Optimization.solve(prob, optimizer, callback = callback, maxiters = 500)

# Train for more iterations with a smaller learning rate
optimizer = OptimizationOptimisers.Adam(0.001)
res = Optimization.solve(prob, res.u, callback = callback, maxiters = 1000)


# Get the trained parameters and the solution function
phi, params = res.u

# Helper functions to get u and v from the NN solution
# phi[1] is the solution for u(x,y)
# phi[2] is the solution for v(x,y)
sol_u(x, y) = phi[1]([x, y], params)[1]
sol_v(x, y) = phi[2]([x, y], params)[1]

# --- Plotting ---
using Plots

# Create a grid of points
xs = 0.0:L/20:L
ys = -h/2:h/10:h/2

# Evaluate solution on the grid
u_sol = [sol_u(x, y) for x in xs, y in ys]
v_sol = [sol_v(x, y) for x in xs, y in ys]

# Plot the vertical displacement v
heatmap(xs, ys, v_sol', xlabel="x (Length)", ylabel="y (Height)",
        title="Vertical Displacement (v)", aspect_ratio=:equal)