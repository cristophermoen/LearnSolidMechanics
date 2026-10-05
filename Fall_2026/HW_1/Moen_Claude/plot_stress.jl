# plot_stress.jl
#
# Stress components on the deformed mid-plane section (y = H/2) for the three
# boundary-condition cases.  One row per case, one column per component.
#
# σxy and σyz are identically zero on y = H/2 because that plane is a plane of
# symmetry (same reason ε_xy and v vanish there, see plot_strain.jl); they are
# shown anyway for completeness, alongside σxz, the in-plane shear produced by
# the barrelling.
#
# Run:  julia --project=. plot_stress.jl

include("cube_compression.jl")
include("section_plot.jl")

results = [solve_cube(c; verbose = false) for (c, _) in CASES]

# stresses are stored in MPa; show colourbars in kPa (σ0 itself is ~1.77 kPa)
comp(a, b) = r -> [1e3 * s[a, b] for s in r.σ_nodes]
comps = [("σ_zz", comp(3, 3), :viridis, false),
         ("σ_xx", comp(1, 1), :viridis, false),
         ("σ_yy", comp(2, 2), :viridis, false),
         ("σ_xy", comp(1, 2), :RdBu,    true),
         ("σ_xz", comp(1, 3), :RdBu,    true),
         ("σ_yz", comp(2, 3), :RdBu,    true)]

title = @sprintf("Stress on section y = H/2   (TPR: H = %.0f mm, E = %.1f MPa, ν = %.2f; 5 oz tuna can: σ0 = %.3f kPa",
                 H, E, ν, 1e3 * σ0) * ")"
plot_component_grid(results, comps, title, "cube_stress.png"; scale = 1e-3, unit = "kPa")
