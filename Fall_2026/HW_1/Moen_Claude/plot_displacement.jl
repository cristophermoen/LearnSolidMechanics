# plot_displacement.jl
#
# Displacement components u, v, w and the total magnitude |u| on the deformed
# mid-plane section (y = H/2) for the three boundary-condition cases.
#
# Note: v is identically zero on y = H/2 (plane of symmetry), so that column is
# blank by physics; v on the side faces y = 0, H equals ∓ the u shown on x = 0, H.
#
# Run:  julia --project=. plot_displacement.jl

include("cube_compression.jl")
include("section_plot.jl")

results = [solve_cube(c; verbose = false) for (c, _) in CASES]

comp(k) = r -> [d[k] for d in r.u_nodes]
comps = [("u",   comp(1),                       :RdBu,    true),
         ("v",   comp(2),                       :RdBu,    true),
         ("w",   comp(3),                       :viridis, false),
         ("|u|", r -> [norm(d) for d in r.u_nodes], :viridis, false)]

title = @sprintf("Displacement on section y = H/2   (TPR: H = %.0f mm, E = %.1f MPa, ν = %.2f; 5 oz tuna can: σ0 = %.3f kPa",
                 H, E, ν, 1e3 * σ0) * ")"
plot_component_grid(results, comps, title, "cube_displacement.png"; scale = 1e-3, unit = "μm")
