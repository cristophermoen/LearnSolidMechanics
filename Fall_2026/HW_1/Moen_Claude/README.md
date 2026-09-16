# Elastic cube under uniform top compression (Ferrite.jl)

28 mm cube of elastic thermoplastic rubber (TPR, Shore A ≈ 60-70: E = 4 MPa, ν = 0.49)
with a 5 oz can of tuna sitting on top.
The can is wider than the cube, so its weight W = 1.39 N is spread uniformly over the
top face as σ0 = W / H² = 1.77 kPa. Bottom face restrained. Linear hexahedral elements, 14 × 14 × 14 mesh.

## Files

| File | Purpose |
|---|---|
| `cube_compression.jl` | Assembles and solves the 3D elasticity problem, prints a summary, writes `.vtu` files |
| `plot_section.jl` | Plots the deformed mid-plane section (y = H/2) coloured by σzz, saves `cube_section.png` |
| `cube_frictionless.vtu`, `cube_clamped.vtu`, `cube_rough_platens.vtu` | ParaView output (displacement `u`, stress `sigma`). Use *Warp By Vector* on `u` |
| `section_plot.jl` | Shared helpers for the row-per-case, column-per-component section plots |
| `plot_strain.jl` | Plots ε_zz, ε_xx, ε_xy, ε_xz on the same section for all three cases, saves `cube_strain.png` |
| `plot_displacement.jl` | Plots u, v, w and \|u\| on the same section, saves `cube_displacement.png` |
| `element_figure.jl` | Draws the 27-node reference element: node numbering, (ξ, η, ζ) axes, shape functions; saves `element_27node.png` |
| `cube_section.png` | Side-by-side deformed sections coloured by σzz |
| `cube_strain.png` | Strain components, one row per boundary-condition case |
| `cube_displacement.png` | Displacement components and magnitude, same layout |
| `element_27node.png` | Reference element, Ferrite node numbering, 1D Lagrange polynomials, shape-function slices |

## Run

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. cube_compression.jl
julia --project=. plot_section.jl
julia --project=. plot_strain.jl
julia --project=. plot_displacement.jl
julia --project=. element_figure.jl
```

ε_xy and v are identically zero on the y = H/2 section because that plane is a plane of
symmetry; ε_xz is plotted alongside as the in-plane shear produced by the barrelling, and
v on the side faces y = 0, H mirrors the u shown on x = 0, H.

## Three boundary-condition cases

* **Frictionless** (w = 0 only, rigid-body modes removed by u = 0 on x = H/2 and
  v = 0 on y = H/2): exact solution is homogeneous uniaxial stress.
  FE reproduces Δw = σ0 H / E = 0.0124 mm, width increase ν σ0 H / E = 0.0061 mm
  and σzz = −1.77 kPa everywhere. This verifies the code.
* **Clamped** (u = v = w = 0): the bottom cannot spread, so the sides barrel outward
  over roughly the lower third of the height, the top shortens slightly less
  (0.0113 mm mean, 91 % of the uniaxial value), and σzz concentrates at the clamped
  corners (≈ −6.7 kPa on the 14³ mesh, growing with refinement because the corner is
  singular in linear elasticity).
* **Rough platens** (bottom clamped, top pressed by a rigid rough platen: u = v = 0 and
  uniform w on the top face): the problem is solved with a trial platen displacement and
  then scaled so the platen force equals W = 1.39 N. Both faces are held, so the
  cube barrels symmetrically, top shortening drops to 0.0101 mm (81 % of uniaxial), and
  σzz concentrates at all eight corners.

## Note on ν = 0.49 and element choice

Rubber is nearly incompressible. ν is set to 0.49 rather than 0.5 because a pure
displacement formulation locks as ν → 0.5. Fully integrated *linear* hexahedra are
already unreliable at 0.49: for the clamped case the mean top shortening creeps from
79 % to 89 % of the uniaxial value between 6³ and 20³ meshes, and σzz oscillates
node-to-node at the clamped face (−3.9, −0.8, −1.8 kPa at z = 0, 2, 4 mm on the 14³
mesh) with a spurious tensile band. *Quadratic* hexahedra (`order = 2`) remove the
oscillation and converge: 90.3 % on 8³ and 90.7 % on 14³. The scripts therefore use
`order = 2`. The frictionless case is exact for either order because its strain field
is uniform.

The plots amplify displacements so the largest one draws as 15 % of H; the factor is
printed in each figure title. Edit the constants at the top of `cube_compression.jl`
(`H`, `E`, `ν`, `m_can`, `nel`, `order`) to change the problem.
