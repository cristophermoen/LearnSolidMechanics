# plot_strain.jl
#
# Strain components on the deformed mid-plane section (y = H/2) for the three
# boundary-condition cases.  One row per case, one column per component.
#
# Note: ε_xy is identically zero on y = H/2 because that plane is a plane of
# symmetry (v is odd in y - H/2, u is even), so ε_xz is shown as well: it is the
# in-plane shear that accompanies the barrelling.
#
# Run:  julia --project=. plot_strain.jl

include("cube_compression.jl")
include("section_plot.jl")

results = [solve_cube(c; verbose = false) for (c, _) in CASES]

# common power-of-ten scale for the strain colourbars
ε_max = maximum(r -> maximum(e -> maximum(abs, e), r.ε_nodes), results)
expo  = floor(Int, log10(ε_max))

comp(a, b) = r -> [e[a, b] for e in r.ε_nodes]
comps = [("ε_zz", comp(3, 3), :viridis, false),
         ("ε_xx", comp(1, 1), :viridis, false),
         ("ε_xy", comp(1, 2), :RdBu,    true),
         ("ε_xz", comp(1, 3), :RdBu,    true)]

title = @sprintf("Strain on section y = H/2   (TPR: H = %.0f mm, E = %.1f MPa, ν = %.2f; 5 oz tuna can: σ0 = %.3f kPa",
                 H, E, ν, 1e3 * σ0) * ")"
plot_component_grid(results, comps, title, "cube_strain.png"; scale = 10.0^expo, unit = pow10_label(expo))
