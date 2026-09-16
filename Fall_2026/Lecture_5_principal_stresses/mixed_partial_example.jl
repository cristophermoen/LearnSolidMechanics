# Run from this folder with:  julia --project=. mixed_partial_example.jl
#
# Mixed partials don't automatically commute. The standard counterexample is
#
#            ⎧ xy(x² - y²)/(x² + y²)    (x, y) ≠ (0, 0)
#   f(x,y) = ⎨
#            ⎩ 0                        (x, y) = (0, 0)
#
# f is continuous with continuous first partials everywhere, but at the origin 
# ∂²f/∂x∂y = -1 while ∂²f/∂y∂x = +1 because the second partials are not continuous there.

using WGLMakie

f(x, y) = (x == 0 && y == 0) ? 0.0 : x * y * (x^2 - y^2) / (x^2 + y^2)

#first partials away from the origin, from the quotient rule 
fx(x, y) = (x == 0 && y == 0) ? 0.0 : y * (x^4 + 4x^2 * y^2 - y^4) / (x^2 + y^2)^2
fy(x, y) = (x == 0 && y == 0) ? 0.0 : x * (x^4 - 4x^2 * y^2 - y^4) / (x^2 + y^2)^2

#second partials away from the origin, from the quotient rule 
#f_xy = f_yx everywhere except the origin, where the mixed partial does not exist as a continuous function
fxx(x, y) = (x == 0 && y == 0) ? 0.0 : -4x * y^3 * (x^2 - 3y^2) / (x^2 + y^2)^3
fyy(x, y) = (x == 0 && y == 0) ? 0.0 : -4x^3 * y * (3x^2 - y^2) / (x^2 + y^2)^3
fxy(x, y) = (x == 0 && y == 0) ? NaN : (x^6 + 9x^4 * y^2 - 9x^2 * y^4 - y^6) / (x^2 + y^2)^3

#check the closed-form second partials against central finite differences at a test point 
x0, y0, h = 0.3, -0.7, 1e-5
@assert isapprox(fxx(x0, y0), (fx(x0 + h, y0) - fx(x0 - h, y0)) / 2h, atol = 1e-5)
@assert isapprox(fyy(x0, y0), (fy(x0, y0 + h) - fy(x0, y0 - h)) / 2h, atol = 1e-5)
@assert isapprox(fxy(x0, y0), (fx(x0, y0 + h) - fx(x0, y0 - h)) / 2h, atol = 1e-5)
@assert isapprox(fxy(x0, y0), (fy(x0 + h, y0) - fy(x0 - h, y0)) / 2h, atol = 1e-5)

#mixed partials at the origin from the limit definition 
h = 1e-6
fxy_0 = (fx(0, h) - fx(0, 0)) / h     # ∂/∂y of f_x at (0,0), expect -1
fyx_0 = (fy(h, 0) - fy(0, 0)) / h     # ∂/∂x of f_y at (0,0), expect +1
println("∂²f/∂y∂x at origin = ", fxy_0)
println("∂²f/∂x∂y at origin = ", fyx_0)


############################################################################
# plot
############################################################################

xs = range(-1, 1, length = 201)
ys = range(-1, 1, length = 201)
z = [f(x, y) for x in xs, y in ys]

fig = Figure(size = (1800, 1200))

#surface of f(x,y)
ax1 = Axis3(fig[1, 1], 
            xlabel = "x", ylabel = "y", zlabel = "f(x,y)",
            title = "f(x,y) = xy(x² - y²)/(x² + y²)",
            azimuth = 0.25π, elevation = 0.2π)
surface!(ax1, xs, ys, z, colormap = :viridis)

#f is zero along both axes and along the diagonals y = ±x 
lines!(ax1, xs, zeros(length(xs)), zeros(length(xs)), color = :red, linewidth = 3)
lines!(ax1, zeros(length(ys)), ys, zeros(length(ys)), color = :red, linewidth = 3)

#first partials along the axes, whose slopes at the origin are the mixed partials 
ax2 = Axis(fig[1, 2:3], 
           xlabel = "distance along axis", ylabel = "first partial derivative",
           title = "f_x(0, y) = -y  and  f_y(x, 0) = x")
lines!(ax2, ys, fx.(0.0, ys), linewidth = 3, label = "f_x(0, y), slope ∂²f/∂y∂x = -1")
lines!(ax2, xs, fy.(xs, 0.0), linewidth = 3, label = "f_y(x, 0), slope ∂²f/∂x∂y = +1")
axislegend(ax2, position = :lt)

#second partial derivatives 
#these are homogeneous of degree zero, so they are bounded but have no limit at the origin, 
#the value there depends on the direction of approach 
#the grid has an even number of points so the origin itself is not sampled 
xs2 = range(-1, 1, length = 200)
ys2 = range(-1, 1, length = 200)
second_partials = [(fxx, "∂²f/∂x²"), (fyy, "∂²f/∂y²"), (fxy, "∂²f/∂x∂y = ∂²f/∂y∂x  (away from origin)")]
for (col, (g, name)) in enumerate(second_partials)
    zg = [g(x, y) for x in xs2, y in ys2]
    axg = Axis3(fig[2, col], 
                xlabel = "x", ylabel = "y", zlabel = "", title = name,
                azimuth = 0.25π, elevation = 0.2π)
    surface!(axg, xs2, ys2, zg, colormap = :balance, colorrange = (-3, 3))
end

#the mixed partial approaches +1 along the x axis and -1 along the y axis 
axm = fig.content[end]
lines!(axm, xs2, zeros(200), fill(1.0, 200), color = :black, linewidth = 3, linestyle = :dash)
lines!(axm, zeros(200), ys2, fill(-1.0, 200), color = :black, linewidth = 3, linestyle = :dash)

Colorbar(fig[2, 4], colormap = :balance, colorrange = (-3, 3), label = "second partial derivative")

display(fig)

#WGLMakie handles mouse rotation of Axis3 on the Julia side, so if this file is run as a 
#script (julia --project=. file.jl) Julia would exit right after display and the browser 
#page would freeze. Keep the process alive until the user presses Enter. 
#Not needed when running in the REPL or VS Code, where Julia stays open.
if !isinteractive()
    println("figure is live in the browser, press Enter here to exit")
    readline()
end
