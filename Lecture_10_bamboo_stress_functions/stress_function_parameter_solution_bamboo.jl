
using CairoMakie, ForwardDiff, StaticArrays, NonlinearSolve

#Airy stress function for axisymmetric case
ϕ(r, a0, a1, a2, a3) = a0 + a1 * log(r) + a2 * r^2 + a3 * r^2 * log(r)

#Radial stress for axisymmetric case 
#σ_r = 1/r * ∂ϕ/∂r + ∂²ϕ/∂θ²
σ_r(r, dϕ_dr) = 1/r * dϕ_dr


σ_r_inner = 0.0
σ_r_outer = 100.0  #N/mm^2
r_inner = 18.0 #mm
r_outer = 25.0 #mm



p = [σ_r_inner, σ_r_outer, r_inner, r_outer]

function calculate_radial_stress_residual(u, p)

    σ_r_inner, σ_r_outer, r_inner, r_outer = p

    a0 = u[1]
    a1 = u[2]
    a2 = u[3]
    a3 = u[4]

    dϕ_dr_r_inner = ForwardDiff.derivative(r -> ϕ(r, a0, a1, a2, a3), r_inner)
    dϕ_dr_r_outer = ForwardDiff.derivative(r -> ϕ(r, a0, a1, a2, a3), r_outer)

    residual_inner = abs(σ_r_inner - σ_r(r_inner, dϕ_dr_r_inner))
    residual_outer = abs(σ_r_outer - σ_r(r_outer, dϕ_dr_r_outer))
    
    residual = residual_inner + residual_outer

    return [residual]

end


u0 = [1.0, 1.0, 1.0, 1.0]
prob = NonlinearProblem(calculate_radial_stress_residual, u0, p)
solver = solve(prob)


calculate_radial_stress_residual(solver.u, p)


function calculate_radial_stress(a0, a1, a2, a3, r_val)

    dϕ_dr = ForwardDiff.derivative(r -> ϕ(r, a0, a1, a2, a3), r_val)

    σ = σ_r(r_val, dϕ_dr)

    return σ

end


r_val = 25.0
a0, a1, a2, a3 = solver.u 

σ_r_outer_confirm = calculate_radial_stress(a0, a1, a2, a3, r_val)


radius_range = range(r_inner, r_outer, 20)

σ_r_range = calculate_radial_stress.(a0, a1, a2, a3, radius_range)

f = Figure()
ax = Axis(f[1, 1], xlabel = "r[mm]",
ylabel = "σr [N/mm^2]")
scatterlines!(radius_range, σ_r_range)
xlims!(ax, 0.0, 25.0)
ylims!(ax, 0.0, 200)
f