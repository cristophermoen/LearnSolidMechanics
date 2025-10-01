
using CairoMakie, ForwardDiff

#Airy stress function 
ϕ(r, a0, a1, a2, a3) = a0 + a1 * log(r) + a2 * r^2 + a3 * r^2 * log(r)

#control parameters, probably need to solve for these using boundary conditions 
a0 = 1.0
a1 = 1.0
a2 = 1.0
a3 = 1.0
rs = 18:0.05:25.0


stress_function = ϕ.(rs, a0, a1, a2, a3)

f = Figure()
ax = Axis(f[1, 1], xlabel = "r",
ylabel = "stress function magnitude")
scatterlines!(rs, stress_function)
xlims!(ax, 0.0, 50.0)
ylims!(ax, 0.0, 3000)
f

#we want to calculate σ_r using the Airy stress function 
#σ_r = 1/r * ∂ϕ/∂r + ∂²ϕ/∂θ²

#∂ϕ/∂r
r_val = 18.0
dϕ_dr_val = ForwardDiff.derivative(r -> ϕ(r, a0, a1, a2, a3), r_val)

σ_r(r, dϕ_dr) = 1/r * dϕ_dr

σ_r(r_val, dϕ_dr_val)



f(u, p) = u .* u .- 2
u0 = @SVector[1.0, 1.0]
prob = NonlinearProblem(f, u0)
solver = solve(prob)





# f = Figure(size = (800, 800))

# ax = PolarAxis(f[1, 1], title = "Airy stress function")
# θs = range(0, 2pi, 37)
# stress_field = [ϕ.(rs, a0, a1, a2, a3) for θ in θs, r in rs]
# p = surface!(ax, 0..2pi, 0.7..1, zeros(size(stress_field)), color = stress_field, shading = NoShading, colormap = :coolwarm)
# f
