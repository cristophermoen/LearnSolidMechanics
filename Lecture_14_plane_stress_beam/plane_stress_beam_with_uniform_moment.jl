using ModelingToolkit, MethodOfLines, DomainSets, NonlinearSolve

@parameters x y
@variables u(..)
Dx = Differential(x)
Dy = Differential(y)

eq = Dx(u(x, y)) + (3 * M)/(2*E*c^3) * y ~ 0

bcs = [u(0, y) ~ 0.0,
       u(L, y) ~ 0.0,
       u(x, 0) ~ x * y,
       u(x, 1) ~ x * y]


# Space and time domains
domains = [x ∈ Interval(0.0, 1.0),
           y ∈ Interval(0.0, 1.0)]

@named pdesys = PDESystem([eq], bcs, domains, [x, y], [u(x, y)])

dx = 0.1
dy = 0.1

# Note that we pass in `nothing` for the time variable `t` here since we
# are creating a stationary problem without a dependence on time, only space.
discretization = MOLFiniteDifference([x => dx, y => dy], nothing, approx_order=2)

prob = discretize(pdesys, discretization)
sol = NonlinearSolve.solve(prob, NewtonRaphson())

u_sol = sol[u(x, y)]

using Plots

heatmap(sol[x], sol[y], u_sol, xlabel="x values", ylabel="y values",
        title="Steady State Heat Equation")