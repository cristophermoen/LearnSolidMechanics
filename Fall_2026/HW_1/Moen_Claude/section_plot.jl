# section_plot.jl
#
# Shared plotting helpers: draw a nodal scalar field on the deformed mid-plane
# section (y = H/2) of the cube, one row per boundary-condition case and one
# column per field component.  Included by plot_strain.jl and plot_displacement.jl.

import CairoMakie                      # qualified: Makie also exports `Vec`
using CairoMakie: Figure, Axis, Label, Colorbar, DataAspect, Point2f, TriangleFace,
                  mesh!, linesegments!, lines!, save, NoShading, hidexdecorations!, hideydecorations!
using Statistics: quantile
using Printf

# displacement amplification for plotting: the max displacement draws as 15% of H
amp_for(res) = 0.15 * H / maximum(norm, res.u_nodes)

# nodal scalar field `vals[node]` -> (deformed points, triangle faces, colours, element edges)
function section_data(res, vals, AMP)
    grid, u_nodes = res.grid, res.u_nodes
    nodes = getnodes(grid)
    ids  = [i for (i, n) in enumerate(nodes) if abs(n.x[2] - H / 2) < 1e-8]
    pos  = Dict(i => k for (k, i) in enumerate(ids))
    Xd   = [Point2f(nodes[i].x[1] + AMP * u_nodes[i][1], nodes[i].x[3] + AMP * u_nodes[i][3]) for i in ids]
    c    = [vals[i] for i in ids]
    faces = TriangleFace{Int}[]
    edges = Point2f[]
    for cell in getcells(grid)
        on = [n for n in cell.nodes if haskey(pos, n)]
        length(on) == 4 || continue
        q = sort(on, by = n -> (round(nodes[n].x[3]; digits = 6), round(nodes[n].x[1]; digits = 6)))
        a, b, cc, d = pos[q[1]], pos[q[2]], pos[q[3]], pos[q[4]]   # BL, BR, TL, TR
        push!(Xd, sum(Xd[[a, b, cc, d]]) / 4)                      # centroid -> symmetric colouring
        push!(c,  sum(c[[a, b, cc, d]]) / 4)
        m = length(Xd)
        push!(faces, TriangleFace(a, b, m)); push!(faces, TriangleFace(b, d, m))
        push!(faces, TriangleFace(d, cc, m)); push!(faces, TriangleFace(cc, a, m))
        append!(edges, (Xd[a], Xd[b], Xd[b], Xd[d], Xd[d], Xd[cc], Xd[cc], Xd[a]))
    end
    return Xd, faces, c, edges
end

function draw_field!(ax, res, vals, crange, cmap, AMP)
    Xd, faces, c, edges = section_data(res, vals, AMP)
    mesh!(ax, Xd, faces; color = c, colormap = cmap, colorrange = crange, shading = NoShading)
    linesegments!(ax, edges; color = (:black, 0.3), linewidth = 0.4)
    lines!(ax, Point2f[(0, 0), (H, 0), (H, H), (0, H), (0, 0)]; color = :gray40, linestyle = :dash, linewidth = 1.2)
end

# unicode superscript for "×10^expo"
function pow10_label(expo)
    sup = Dict('-' => '⁻', '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
               '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹')
    return "×10" * join(sup[c] for c in string(expo))
end

const CASES = [(:frictionless,  "Frictionless bottom\nfree top"),
               (:clamped,       "Clamped bottom\ntraction on top"),
               (:rough_platens, "Clamped bottom\nrough platen on top")]

# comps: vector of (name, getter(res) -> Vector{Float64} per node, colormap, symmetric?)
# scale: divide values by `scale` for the colourbar ticks; `unit` is the tick unit label
function plot_component_grid(results, comps, title, filename; scale = 1.0, unit = "")
    AMP = amp_for(results[1])
    fig = Figure(size = (400 * length(comps) + 120, 380 * length(CASES) + 60), fontsize = 13)
    Label(fig[0, 1:length(comps)], title * @sprintf("; displacements ×%.3g", AMP),
          fontsize = 15, tellwidth = false)
    for (j, (name, getter, cmap, symm)) in enumerate(comps)
        fields = [getter(r) for r in results]
        allvals = vcat(fields...)
        lo, hi = quantile(allvals, 0.02), quantile(allvals, 0.98)   # clip corner singularities
        if symm
            m = max(abs(lo), abs(hi), 1e-12); lo, hi = -m, m
        elseif hi - lo < 1e-12 * max(abs(hi), 1.0)
            pad = 0.05 * max(abs(hi), 1e-12); lo -= pad; hi += pad
        end
        for (i, r) in enumerate(results)
            ax = Axis(fig[i, j]; aspect = DataAspect(),
                      title = i == 1 ? name : "",
                      ylabel = j == 1 ? CASES[i][2] : "", ylabelsize = 13,
                      xlabel = i == length(CASES) ? "x [mm]" : "")
            i < length(CASES) && hidexdecorations!(ax; grid = false)
            j > 1 && hideydecorations!(ax; grid = false)
            draw_field!(ax, r, fields[i], (lo, hi), cmap, AMP)
        end
        Colorbar(fig[length(CASES) + 1, j]; limits = (lo, hi), colormap = cmap, vertical = false,
                 label = "$name  [$unit]", flipaxis = false,
                 tickformat = v -> [string(round(x / scale; digits = 2)) for x in v])
    end
    CairoMakie.colgap!(fig.layout, 28)          # keep neighbouring colourbar ticks apart
    save(filename, fig; px_per_unit = 2)
    println("wrote $filename")
end
