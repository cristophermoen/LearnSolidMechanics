# Thin plate with a circular hole under uniaxial tension — Airy stress function solution

Reference solution for the **Neural Operator + FEA** in-class project (Fall 2026).

Plate 300 mm × 50 mm with a central hole of diameter 25 mm, pulled at its ends by a uniform
tension `T`; linear elastic `E, ν`; plane stress. The closed-form solution is Kirsch's (1898)
infinite-plate solution, Sadd *Elasticity* Example 8-7, built here from the Michell terms of the
polar Airy stress function with `Symbolics.jl`.

`plate_with_hole_airy.jl` is a [Pluto](https://plutojl.org) notebook that

1. explains which parts of the figure can be solved in closed form (the hole, the ends) and which
   cannot (the long edges, only two hole radii away);
2. writes the plane-stress governing equations in polar coordinates and the boundary conditions,
   with the far-field tension derived from `φ₀ = ½Ty²`;
3. takes the trial function `φ = a₁ ln r + a₂r² + a₃r² ln r + (b₁r² + b₂r⁴ + b₃/r² + b₄) cos 2θ`
   (Sadd 8.4.12), imposes boundedness at infinity, the far field and the traction-free hole, and
   solves for the coefficients exactly (Sadd 8.4.14);
4. recovers the Kirsch stresses (Sadd 8.4.15), the hoop stress `T(1 − 2cos 2θ)` on the hole,
   `K_t = 3`, the net-section profile `σx(0,y) = T(1 + a²/2y² + 3a⁴/2y⁴)`, and the Cartesian
   components, and checks every boundary condition symbolically;
5. integrates the strains to the displacements, verifies them against the Kolosov-constant form in
   the literature, and finds the deformed hole: an ellipse with axis changes `+6Ta/E` and `−2Ta/E`,
   independent of `ν`;
6. evaluates everything for `T = 100 MPa`, `a = 12.5 mm`, `W = 50 mm`, `L = 300 mm`, `E = 200 GPa`,
   `ν = 0.3`, plots the fields on the strip and the deformed shape, plots the residual tractions the
   infinite-plate solution leaves on the strip edges, shows that the Kirsch field carries only 84 %
   of the load through the strip's net section, and quotes the finite-width (Howland / Peterson)
   stress concentration factor, `K_tg ≈ 4.31` against Kirsch's 3;
7. provides `kirsch(x, y)` to sample all fields on any grid, plus a list of checks to make against
   the finite element model.

## Running

```julia
using Pluto; Pluto.run()
```

then open `plate_with_hole_airy.jl`. Pluto installs `Symbolics`, `CairoMakie` and `Latexify` into
the notebook's own environment on first run (a few minutes).

## Output

| `stress_lines.png` | hoop stress on the hole, net-section and load-axis profiles |
|---|---|
| `stress_fields.png` | `σx, σy, τxy` on the strip and near the hole |
| `edge_residuals.png` | tractions the infinite-plate solution leaves on the strip edges and ends |
| `deformed_shape.png` | deformed grid and the elliptical hole, displacements × 200 |
