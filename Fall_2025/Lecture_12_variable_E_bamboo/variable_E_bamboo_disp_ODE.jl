using DifferentialEquations, NonlinearSolve, Polynomials, ForwardDiff


Eo = 10000.0 #N/mm^2
ν = 0.20
r_inner = 18.0 #mm 
r_outer = 25.0 #mm 

σ_r_inner = 0.0 #N/mm^2 
σ_r_outer = -20.0 #N/mm^2 

a = r_inner #mm, inner radius 
n = 1 #linear variation 

E(r, a, n) = Eo * (r / a) ^ n # Elastic modulus varying through bamboo thickness 

μ(r, a, n) = E(r, a, n) / (2 * (1 + ν)) #shear modulus

λ(r, a, n) = (ν * E(r, a, n)) / (2 * (1 + ν)) #Lame's constant   

function navier_equation(urr, ur, u, p, r)

    #urr is the second derivative of u with respect to r 
    #ur is the first deriviative of u with respect to r 
    #u is u(r) radial deformation through the bamboo thickness 
	
    n, ν = p

    #Eq. 14.2.4 in Sadd 
    urr[1] = -(n+1)/r * ur[1] - (n * ν - 1) * u[1] / r^2

end


p = [n, ν, a, r_inner, r_outer, σ_r_inner, σ_r_outer]



function calculate_radial_stress_residual(u, p)

    n, ν, a, r_inner, r_outer, σ_r_inner, σ_r_outer = p

# ur_0 = [0.005]
# u_0 = [-1.17]


    ur_0 = [u[1]] 
    u_0 = [u[2]]

    prob = SecondOrderODEProblem(navier_equation, ur_0, u_0, (r_inner, r_outer), p);

    solution = solve(prob, DPRKN6());

    r_for_solution = solution.t 
    u_for_solution = [solution.u[i][2] for i in eachindex(solution.u)]


    u_fit = fit(r_for_solution, u_for_solution)

    ur_fit = derivative(u_fit)

    ϵ_r_inner = ur_fit(r_inner)
    ϵ_θ_inner = 1/r_inner * u_fit(r_inner)
    σ_r_inner_solution = λ(r_inner, a, n) * (ϵ_r_inner + ϵ_θ_inner) + 2 * μ(r_inner, a, n) * ϵ_r_inner


    ϵ_r_outer = ur_fit(r_inner)
    ϵ_θ_outer = 1/r_outer * u_fit(r_outer)
    σ_r_outer_solution = λ(r_outer, a, n) * (ϵ_r_outer +  ϵ_θ_outer) + 2 * μ(r_outer, a, n) * ϵ_r_outer


    residual_inner = abs(σ_r_inner - σ_r_inner_solution)
    residual_outer = abs(σ_r_outer - σ_r_outer_solution)

    residual = residual_inner + residual_outer 

    return [residual]

end


u0 = [0.1, 0.1]
prob = NonlinearProblem(calculate_radial_stress_residual, u0, p)
solver = solve(prob)

### calculate radial and tangential strains


    u_fit = fit(r_for_solution, u_for_solution)

    ur_fit = derivative(u_fit)

    ϵ_r_inner = ur_fit(r_inner)
    ϵ_θ_inner = 1/r_inner * u_fit(r_inner)
 


###calculate radial and tangenetial stresses 

   σ_r_inner_solution = λ(r_inner, a, n) * (ϵ_r_inner + ϵ_θ_inner) + 2 * μ(r_inner, a, n) * ϵ_r_inner


