# element_figure.jl
#
# The 27-node quadratic Lagrange hexahedron used in cube_compression.jl, drawn
# in its reference (unit) configuration with Ferrite's node numbering, the
# reference coordinate system (ξ, η, ζ) ∈ [-1, 1]³, and its shape functions.
#
# Run:  julia --project=. element_figure.jl

using Ferrite, LinearAlgebra, Printf
import CairoMakie
using CairoMakie: Figure, Axis, Axis3, Label, Colorbar, Legend, MarkerElement, DataAspect,
                  GridLayout, Point3f, Point2f, scatter!, lines!, text!, heatmap!, hidedecorations!,
                  hidespines!, save, Relative, colsize!, rowsize!, rowgap!, colgap!

ip = Lagrange{RefHexahedron,2}()
X  = Ferrite.reference_coordinates(ip)              # 27 node positions in (ξ, η, ζ)
nnodes = length(X)

# --- node classification (Ferrite ordering) ---------------------------------
kind(i) = i <= 8 ? :corner : i <= 20 ? :edge : i <= 26 ? :face : :interior
kcolor  = Dict(:corner => :firebrick, :edge => :royalblue3, :face => :seagreen, :interior => :darkorange)
klabel  = Dict(:corner => "corner  (1–8)", :edge => "edge midpoint  (9–20)",
               :face => "face centre  (21–26)", :interior => "centroid  (27)")

# --- 1D quadratic Lagrange polynomials on [-1, 1] ----------------------------
L(t, a) = a == -1 ? t * (t - 1) / 2 : a == 0 ? 1 - t^2 : t * (t + 1) / 2
N(i, ξ, η, ζ) = L(ξ, Int(X[i][1])) * L(η, Int(X[i][2])) * L(ζ, Int(X[i][3]))

# ============================================================================
fig = Figure(size = (1750, 1300), fontsize = 14)
Label(fig[0, 1], "27-node quadratic Lagrange hexahedron (Ferrite Lagrange{RefHexahedron,2}), reference configuration",
      fontsize = 20, font = :bold, tellwidth = false)
top    = fig[1, 1] = GridLayout()         # [3D element | 1D polys + formula | table]
right  = top[1, 2] = GridLayout()
bottom = fig[2, 1] = GridLayout()         # shape-function slices

# --- (A) 3D reference element -------------------------------------------------
ax = Axis3(top[1, 1]; aspect = :data, azimuth = 1.33π, elevation = 0.10π,
           title = "Reference element  (ξ, η, ζ) ∈ [-1, 1]³   (nodes numbered as in Ferrite)",
           titlesize = 16, viewmode = :fit,
           limits = ((-1.7, 1.7), (-1.7, 1.7), (-1.7, 1.7)))
hidedecorations!(ax); hidespines!(ax)

# 12 edges: (corner a, corner b) pairs from the edge dof tuples
for e in Ferrite.edgedof_indices(ip)
    a, b = X[e[1]], X[e[2]]
    lines!(ax, [Point3f(a...), Point3f(b...)]; color = :gray35, linewidth = 2)
end
# dotted mid-planes ζ = 0 and the face-centre plane through the centroid, to place interior nodes
for pts in ([(-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0),(-1,-1,0)],
            [(0,-1,-1),(0,1,-1),(0,1,1),(0,-1,1),(0,-1,-1)],
            [(-1,0,-1),(1,0,-1),(1,0,1),(-1,0,1),(-1,0,-1)])
    lines!(ax, [Point3f(p...) for p in pts]; color = (:gray60, 0.6), linestyle = :dot, linewidth = 1)
end

# nodes
for k in (:corner, :edge, :face, :interior)
    ids = [i for i in 1:nnodes if kind(i) == k]
    scatter!(ax, [Point3f(X[i]...) for i in ids]; color = kcolor[k], markersize = 22,
             strokecolor = :black, strokewidth = 1)
end
# node numbers, offset outward from the centroid so they don't sit on the marker
for i in 1:nnodes
    p  = X[i]
    d  = norm(p) > 0 ? p / norm(p) : Vec(0.0, 0.0, 1.0)
    q  = p + 0.17 * d
    text!(ax, Point3f(q...); text = string(i), fontsize = 15, font = :bold, align = (:center, :center))
end

# reference axes ξ, η, ζ drawn from the centroid
# labels are offset off the axis line so they clear the face-centre nodes 23, 24, 26
for (v, name, off) in ((Vec(1.0,0.0,0.0), "ξ", Vec(0.0, 0.0, 0.16)),
                       (Vec(0.0,1.0,0.0), "η", Vec(0.0, 0.0, 0.16)),
                       (Vec(0.0,0.0,1.0), "ζ", Vec(0.14, 0.0, 0.0)))
    lines!(ax, [Point3f(0,0,0), Point3f((1.5 .* v)...)]; color = :black, linewidth = 2.5)
    text!(ax, Point3f((1.5 .* v .+ off)...); text = name, fontsize = 22, font = :bold_italic, align = (:center, :center))
end

Legend(top[1, 1], [MarkerElement(color = kcolor[k], marker = :circle, markersize = 18, strokecolor = :black, strokewidth = 1)
                     for k in (:corner, :edge, :face, :interior)],
       [klabel[k] for k in (:corner, :edge, :face, :interior)];
       tellheight = false, tellwidth = false, halign = :left, valign = :top, margin = (10, 10, 10, 10),
       framevisible = true, backgroundcolor = (:white, 0.85))

# --- (B) 1D Lagrange polynomials ----------------------------------------------
ax1 = Axis(right[1, 1]; title = "1D quadratic Lagrange polynomials on [-1, 1]", titlesize = 15,
           xlabel = "t", ylabel = "L(t)", xticks = -1:0.5:1)
t = range(-1, 1; length = 201)
lines!(ax1, t, L.(t, -1); linewidth = 3, label = "L₋(t) = t(t − 1)/2   (node at t = −1)")
lines!(ax1, t, L.(t,  0); linewidth = 3, label = "L₀(t) = 1 − t²          (node at t =  0)")
lines!(ax1, t, L.(t, +1); linewidth = 3, label = "L₊(t) = t(t + 1)/2   (node at t = +1)")
scatter!(ax1, [-1, 0, 1], [1, 1, 1]; color = :black, markersize = 9)
scatter!(ax1, [-1, 0, 1, -1, 0, 1], [0, 0, 0, 0, 0, 0]; color = :black, markersize = 9)
Legend(right[2, 1], ax1; orientation = :horizontal, framevisible = false, labelsize = 12, tellheight = true)

# --- (C) formula + properties ---------------------------------------------------
formula = """
Shape function of node i with reference coordinates (ξᵢ, ηᵢ, ζᵢ), each ∈ {−1, 0, +1}:

        Nᵢ(ξ, η, ζ) = L_ξᵢ(ξ) · L_ηᵢ(η) · L_ζᵢ(ζ)          (L₋, L₀, L₊ plotted above)

  •  Nᵢ(ξⱼ, ηⱼ, ζⱼ) = δᵢⱼ    interpolatory: 1 at its own node, 0 at the other 26
  •  Σᵢ Nᵢ = 1               partition of unity: reproduces a rigid translation
  •  triquadratic: complete to degree 2, plus all mixed terms up to ξ²η²ζ²
  •  27 nodes × 3 displacement components = 81 dofs per element

Examples:   N₁  = (1/8) ξ(ξ−1) η(η−1) ζ(ζ−1)     corner        (−1, −1, −1)
            N₉  = (1/4) (1−ξ²) η(η−1) ζ(ζ−1)     edge midpoint ( 0, −1, −1)
            N₂₁ = (1/2) (1−ξ²)(1−η²) ζ(ζ−1)      face centre   ( 0,  0, −1)
            N₂₇ = (1−ξ²)(1−η²)(1−ζ²)             centroid      ( 0,  0,  0)

Displacement in the element:  u(ξ,η,ζ) = Σᵢ Nᵢ(ξ,η,ζ) uᵢ.  The geometry x(ξ,η,ζ)
is mapped with the same Nᵢ (isoparametric); stiffness uses 3×3×3 Gauss points.
"""
Label(right[3, 1]; text = formula, font = "DejaVu Sans Mono", fontsize = 12, justification = :left,
      halign = :left, valign = :top, tellwidth = false, tellheight = true)

# --- (D) node coordinate table ------------------------------------------------
rows = ["  i   ξᵢ  ηᵢ  ζᵢ"]
for i in 1:nnodes
    push!(rows, @sprintf("%3d  %+d  %+d  %+d", i, Int(X[i][1]), Int(X[i][2]), Int(X[i][3])))
end
cols = [join(rows[[1; (2:10)]], "\n"), join(rows[[1; (11:19)]], "\n"), join(rows[[1; (20:28)]], "\n")]
Label(top[1, 3]; text = "Node coordinates\n\n" * cols[1] * "\n\n" * cols[2] * "\n\n" * cols[3],
      font = "DejaVu Sans Mono", fontsize = 13, justification = :left, halign = :left, valign = :top,
      tellwidth = true, tellheight = false)
colsize!(top, 1, Relative(0.46))

# --- (E) shape functions on a slice through their node --------------------------
ξs = range(-1, 1; length = 121)
examples = [(1, "N₁  corner (−1,−1,−1),  slice ζ = −1"),
            (9, "N₉  edge midpoint (0,−1,−1),  slice ζ = −1"),
            (21, "N₂₁  face centre (0,0,−1),  slice ζ = −1"),
            (27, "N₂₇  centroid (0,0,0),  slice ζ = 0")]
hm = nothing
for (j, (i, ttl)) in enumerate(examples)
    ζi = X[i][3]
    Z  = [N(i, ξ, η, ζi) for ξ in ξs, η in ξs]
    axh = Axis(bottom[1, j]; aspect = DataAspect(), title = ttl, titlesize = 13,
               xlabel = "ξ", ylabel = j == 1 ? "η" : "", xticks = -1:1, yticks = -1:1)
    global hm = heatmap!(axh, ξs, ξs, Z; colormap = :RdBu, colorrange = (-1, 1))
    lines!(axh, [-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1]; color = :black, linewidth = 1.5)
    ids = [k for k in 1:nnodes if X[k][3] == ζi]
    scatter!(axh, [X[k][1] for k in ids], [X[k][2] for k in ids]; color = :white, strokecolor = :black,
             strokewidth = 1, markersize = 10)
    scatter!(axh, [X[i][1]], [X[i][2]]; color = kcolor[kind(i)], strokecolor = :black, strokewidth = 1.5, markersize = 16)
end
Colorbar(bottom[1, 5], hm; label = "Nᵢ", ticks = -1:0.5:1)
Label(fig[3, 1], "Bottom row: each shape function equals 1 at its own node (large marker) and 0 at every other node (small markers) lying in the slice; " *
                   "it takes negative values between corner/edge nodes, as quadratic interpolation must.",
      fontsize = 12, tellwidth = false)

rowsize!(fig.layout, 2, Relative(0.30))
save("element_27node.png", fig; px_per_unit = 2)
println("wrote element_27node.png")
