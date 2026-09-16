using ForwardDiff, NonlinearSolve 

μ = 44000.0 #MPa # Ti-6Al-4V
α = 1.0
a = 15.0 / 2
b = 15.0 / 2

ϕ(x, y, a, b, K) = K * (x^2/a^2 + y^2/b^2 - 1)


dϕ2_dx2(K, a) = 2 * K / a^2 
dϕ2_dy2(K, b) = 2 * K / b^2 

K = 1.0
x_range = range(-a, a, 10)
y_range = range(-b, b, 10)

p = (a, b, α, μ)


function calculate_residual(u, p)

    a, b, α, μ = p

    K = u[1]

    residual = dϕ2_dx2(K, a) + dϕ2_dy2(K, b) + 2 * μ * α

    return [residual]

end


# K
u0 = [1.0]
prob = NonlinearProblem(calculate_residual, u0, p)
solver = solve(prob)


#K 