# plot_section.jl
#
# Plot the deformed mid-plane section (y = H/2) of the compressed cube for the
# frictionless and clamped bottom, coloured by σzz, with displacements amplified.
#
# Run:  julia --project=. plot_section.jl

include("cube_compression.jl")
import CairoMakie                      # qualified: Makie also exports `Vec`
using CairoMakie: Figure, Axis, Label, Colorbar, DataAspect, Point2f, TriangleFace, GeometryBasics,
                  mesh!, linesegments!, lines!, save, NoShading
using Statistics: quantile

# displacement amplification for plotting: chosen so the max displacement draws as 15% of H
amp_for(res) = 0.15 * H / maximum(norm, res.u_nodes)
AMP = 1.0   # set after the solves

function section_data(res)
    grid, u_nodes, σ_nodes = res.grid, res.u_nodes, res.σ_nodes
    nodes = getnodes(grid)
    # nodes lying on the plane y = H/2
    ids  = [i for (i, n) in enumerate(nodes) if abs(n.x[2] - H / 2) < 1e-8]
    pos  = Dict(i => k for (k, i) in enumerate(ids))
    X0   = [Point2f(nodes[i].x[1], nodes[i].x[3]) for i in ids]
    Xd   = [Point2f(nodes[i].x[1] + AMP * u_nodes[i][1], nodes[i].x[3] + AMP * u_nodes[i][3]) for i in ids]
    szz  = [σ_nodes[i][3, 3] for i in ids]
    # quads: each hexahedron with 4 nodes on the plane contributes a quad.
    # Split each quad into 4 triangles about its centroid so that the colour
    # interpolation is symmetric (a single diagonal would bias the contours).
    X0 = collect(X0); Xd = collect(Xd); szz = collect(szz)
    faces = TriangleFace{Int}[]
    edges = Point2f[]                       # element edges for the wireframe
    for cell in getcells(grid)
        on = [n for n in cell.nodes if haskey(pos, n)]
        length(on) == 4 || continue
        q = sort(on, by = n -> (round(nodes[n].x[3]; digits = 6), round(nodes[n].x[1]; digits = 6)))
        a, b, c, d = pos[q[1]], pos[q[2]], pos[q[3]], pos[q[4]]   # BL, BR, TL, TR
        push!(X0, sum(X0[[a, b, c, d]]) / 4)
        push!(Xd, sum(Xd[[a, b, c, d]]) / 4)
        push!(szz, sum(szz[[a, b, c, d]]) / 4)
        m = length(Xd)
        push!(faces, TriangleFace(a, b, m)); push!(faces, TriangleFace(b, d, m))
        push!(faces, TriangleFace(d, c, m)); push!(faces, TriangleFace(c, a, m))
        append!(edges, (Xd[a], Xd[b], Xd[b], Xd[d], Xd[d], Xd[c], Xd[c], Xd[a]))
    end
    return X0, Xd, faces, szz, edges
end

function draw!(ax, res, title, crange)
    X0, Xd, faces, szz, edges = section_data(res)
    mesh!(ax, Xd, faces; color = 1e3 .* szz, colormap = :viridis, colorrange = crange, shading = NoShading)
    linesegments!(ax, edges; color = (:black, 0.3), linewidth = 0.5)
    lines!(ax, Point2f[(0, 0), (H, 0), (H, H), (0, H), (0, 0)]; color = :gray40, linestyle = :dash, linewidth = 1.5)
    ax.title = title
end

res_free  = solve_cube(:frictionless;  verbose = false)
res_clamp = solve_cube(:clamped;       verbose = false)
res_plat  = solve_cube(:rough_platens; verbose = false)
AMP = amp_for(res_free)

# colour range in kPa, shared by the three panels (clip the corner singularities)
szz_all = 1e3 .* vcat([[s[3, 3] for s in r.σ_nodes] for r in (res_free, res_clamp, res_plat)]...)
crange = (quantile(szz_all, 0.04), quantile(szz_all, 0.96))   # clip the corner singularities
fig = Figure(size = (1450, 540), fontsize = 14)
Label(fig[0, 1:3],
      @sprintf("TPR cube, H = %.0f mm, E = %.1f MPa, ν = %.2f.  5 oz tuna can on top: W = %.3f N, σ0 = %.3f kPa.  Section y = H/2, displacements ×%.3g",
               H, E, ν, W_can, 1e3 * σ0, AMP),
      fontsize = 15, tellwidth = false)
ax1 = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x [mm]", ylabel = "z [mm]")
ax2 = Axis(fig[1, 2], aspect = DataAspect(), xlabel = "x [mm]")
ax3 = Axis(fig[1, 3], aspect = DataAspect(), xlabel = "x [mm]")
draw!(ax1, res_free,  "Frictionless bottom, free top\nexact uniaxial state",   crange)
draw!(ax2, res_clamp, "Clamped bottom, traction on top\nbarrelling at base",   crange)
draw!(ax3, res_plat,  "Clamped bottom, rigid rough platen on top\nsame force σ0 H²", crange)
Colorbar(fig[1, 4], limits = crange, colormap = :viridis, label = "σzz [kPa]")
save("cube_section.png", fig; px_per_unit = 2)
println("wrote cube_section.png")
