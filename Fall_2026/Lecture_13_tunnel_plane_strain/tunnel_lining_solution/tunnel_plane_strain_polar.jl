### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ d4000000-0000-0000-0000-000000000001
begin
	using Symbolics
	using CairoMakie
	using Latexify      # lets Pluto render Symbolics expressions as LaTeX
end

# ╔═╡ d4000000-0000-0000-0000-000000000002
md"""
# Tunnel Lining Under Uniform Ground Pressure — Plane Strain in Polar Coordinates

A long **semicircular tunnel lining** of inner radius ``R`` and thickness ``t`` (outer radius ``R + t``), made of a linear elastic material ``E, \nu``, is buried so that the surrounding ground presses on its **outer surface** with a uniform radial pressure ``p``. The inner surface is free. The two feet of the arch rest on **supports** at the level of the centre of the arch. The tunnel is long in the direction of its axis, so every cross-section deforms in the same way and the axial strain is zero: **plane strain** (Sadd §7.1).

The geometry is circular, so this is the first problem of the course that calls for the **polar-coordinate** form of the Airy stress function (Sadd §7.6) and the general **Michell solution** (Sadd §8.3). It is the plane-strain twin of Sadd's Example 8-6, the thick-walled cylinder under uniform pressure, with the extra question of what the supports at the feet of the arch do to that solution.

References: Sadd, *Elasticity*, §7.1 (plane strain), §7.2 and Table 7-1 (plane stress ↔ plane strain), §7.6 (polar coordinate formulation), §8.3.1–8.3.2 (Michell solution, axisymmetric case), Example 8-6 and §8.4.1–8.4.2 (thick-walled cylinder and its limiting cases).

**Coordinates.** Origin at the centre of the arch, on the line of the supports; ``x`` horizontal, ``y`` upward; polar angle ``\theta`` measured from the ``+x`` axis, so the right foot is at ``\theta = 0``, the crown at ``\theta = \pi/2`` and the left foot at ``\theta = \pi``. The lining occupies ``R \le r \le R + t``, ``0 \le \theta \le \pi``. The radial displacement ``u_r`` is positive outward and ``u_\theta`` positive counter-clockwise — so at the right foot ``u_\theta`` is the *upward* displacement and at the left foot the *downward* one.

Steps:
1. Write the governing equations in polar coordinates: stresses from ``\phi(r,\theta)``, the biharmonic equation, plane strain Hooke's law, strain–displacement relations.
2. State the boundary conditions, use symmetry, and be explicit about what the supports are taken to do.
3. Pick the terms of the Michell solution that fit the problem and integrate the general displacement field.
4. Enforce the support conditions, then the pressure conditions, symbolically with `Symbolics.jl`.
5. Recover the stresses (including the out-of-plane ``\sigma_z``), the strains and the displacements, and compare with the thin-ring formula of strength of materials.
6. Work a numerical example for a concrete lining with `Makie.jl` plots.
7. Replace the smooth supports by **fixed** ones, see why that problem has no closed-form solution, and solve it in the Saint-Venant sense — which turns out to be the force method of indeterminate arches with exact flexibilities.
"""

# ╔═╡ d4000000-0000-0000-0000-000000000003
md"""
## 1. Governing equations in polar coordinates

With no body forces, an Airy stress function ``\phi(r, \theta)`` generates stresses that satisfy the polar equilibrium equations (Sadd 7.6.3) identically (Sadd 7.6.7):

```math
\sigma_r = \frac{1}{r}\frac{\partial \phi}{\partial r} + \frac{1}{r^2}\frac{\partial^2 \phi}{\partial \theta^2}, \qquad
\sigma_\theta = \frac{\partial^2 \phi}{\partial r^2}, \qquad
\tau_{r\theta} = -\frac{\partial}{\partial r}\left(\frac{1}{r}\frac{\partial \phi}{\partial \theta}\right)
```

Compatibility again reduces to the biharmonic equation, now built from the polar Laplacian (Sadd 7.6.8):

```math
\nabla^4 \phi = \left(\frac{\partial^2}{\partial r^2} + \frac{1}{r}\frac{\partial}{\partial r} + \frac{1}{r^2}\frac{\partial^2}{\partial \theta^2}\right)
\left(\frac{\partial^2 \phi}{\partial r^2} + \frac{1}{r}\frac{\partial \phi}{\partial r} + \frac{1}{r^2}\frac{\partial^2 \phi}{\partial \theta^2}\right) = 0
```

This equation is the same for plane stress and plane strain (Sadd §7.5): the **stresses** of a traction boundary-value problem do not care which one we are in, and do not depend on ``E`` or ``\nu``. The difference appears when we go from stresses to strains. **Plane strain** (``\varepsilon_z = \gamma_{rz} = \gamma_{\theta z} = 0``) Hooke's law, inverted (Sadd 7.1.3 and 7.6.2):

```math
\varepsilon_r = \frac{1+\nu}{E}\left[(1-\nu)\sigma_r - \nu\sigma_\theta\right], \qquad
\varepsilon_\theta = \frac{1+\nu}{E}\left[(1-\nu)\sigma_\theta - \nu\sigma_r\right], \qquad
\gamma_{r\theta} = \frac{2(1+\nu)}{E}\tau_{r\theta}
```

and the out-of-plane normal stress needed to hold ``\varepsilon_z = 0`` is

```math
\sigma_z = \nu\left(\sigma_r + \sigma_\theta\right)
```

Equivalently, take the plane stress formulas of Lectures 11–12 and replace ``E \to E/(1-\nu^2)``, ``\nu \to \nu/(1-\nu)`` (Sadd Table 7-1). Finally, the strain–displacement relations in polar coordinates (Sadd 7.6.1):

```math
\varepsilon_r = \frac{\partial u_r}{\partial r}, \qquad
\varepsilon_\theta = \frac{1}{r}\left(u_r + \frac{\partial u_\theta}{\partial \theta}\right), \qquad
\gamma_{r\theta} = \frac{1}{r}\frac{\partial u_r}{\partial \theta} + \frac{\partial u_\theta}{\partial r} - \frac{u_\theta}{r}
```

Two terms are new compared with Cartesian coordinates. The ``u_r/r`` in the hoop strain says that a ring which simply moves outward gets longer; the ``-u_\theta/r`` in the shear strain is what makes a rigid rotation ``u_\theta = \omega_0 r`` strain-free.
"""

# ╔═╡ d4000000-0000-0000-0000-000000000004
md"""
## 2. Boundary conditions — symmetry, and what the supports are taken to do

**Outer surface** ``r = R + t``, for all ``0 \le \theta \le \pi``: the ground pushes radially inward with intensity ``p`` and (we assume) exerts no shear:

```math
\sigma_r(R + t, \theta) = -p, \qquad \tau_{r\theta}(R + t, \theta) = 0
```

**Inner surface** ``r = R``: traction-free,

```math
\sigma_r(R, \theta) = 0, \qquad \tau_{r\theta}(R, \theta) = 0
```

**The feet** ``\theta = 0`` and ``\theta = \pi``, for all ``R \le r \le R + t``. The sketch shows the feet bearing on supports, and we have to decide what those supports do. The choice that makes the problem solvable in closed form is a **smooth, rigid support**: it stops the foot from moving vertically, lets it slide horizontally, and exerts no friction:

```math
u_\theta(r, 0) = 0, \quad \tau_{r\theta}(r, 0) = 0, \qquad\qquad u_\theta(r, \pi) = 0, \quad \tau_{r\theta}(r, \pi) = 0
```

The vertical reaction is **not** prescribed: it is whatever ``\sigma_\theta(r, 0)`` the solution produces, and we will integrate it afterwards to check that the supports carry the load.

**Symmetry — two of them.** The problem is mirror-symmetric about the vertical axis ``\theta = \pi/2``, so the arch cannot drift sideways as a whole: ``u_\theta(r, \pi/2) = 0``. More useful is a second, hidden symmetry. Imagine the **complete ring** ``0 \le \theta < 2\pi`` under the same uniform external pressure ``p``. It is symmetric about the horizontal diameter, and on a plane of symmetry the displacement normal to the plane (here the vertical one, ``u_\theta``) and the shear traction ``\tau_{r\theta}`` both vanish — *exactly* our support conditions. So the upper half of the full ring **is** the arch on smooth supports, with the lower half playing the part of the supports. And a full ring under uniform pressure is **axisymmetric**: nothing can depend on ``\theta``. That is what lets us keep only the ``\theta``-independent terms of the Michell solution. We will not *assume* axisymmetric displacements, though — the support conditions will enforce it, and in doing so they say something about a term of the general solution that the full ring would never let us see.

### What is relaxed, and what it costs

* **Smooth supports.** A real invert, footing or rock bed resists the small inward slide of the feet by friction, or is cast integral with the lining. The difference between the real and the idealised support is a self-equilibrated horizontal traction on the foot; Saint-Venant's principle (Sadd §5.6) confines its effect to within about one lining thickness ``t`` of the feet. The hoop stress in the haunches and crown is unaffected by that local difference — but a support that *holds* the feet also changes the global bending of the arch, and Section 11 works that case out.
* **Uniform radial pressure.** The ground pressure on a real tunnel is larger vertically than horizontally (the lateral earth-pressure coefficient is usually below one) and is not purely radial. The non-axisymmetric part needs the ``\cos 2\theta`` terms of the Michell solution — the same terms that solve the stress concentration at a hole (Sadd §8.4.3). Here we solve for the axisymmetric part only.
* **Plane strain.** Setting ``\varepsilon_z = 0`` requires a longitudinal stress ``\sigma_z = \nu(\sigma_r + \sigma_\theta)`` along the tunnel. Far from the portals the surrounding ground supplies it; near a free end it is missing and a Saint-Venant end effect appears (Sadd §7.1 describes the corrective solution).
"""

# ╔═╡ d4000000-0000-0000-0000-000000000005
md"""
## 3. Symbolic setup

Independent variables, problem parameters and the differential operators in ``r`` and ``\theta``.
"""

# ╔═╡ d4000000-0000-0000-0000-000000000006
begin
	@variables r θ
	@variables p R t E ν           # ground pressure, inner radius, lining thickness, Young's modulus, Poisson's ratio
	Dr = Differential(r)
	Dθ = Differential(θ)
end

# ╔═╡ d4000000-0000-0000-0000-000000000007
"""
    laplacian_polar(f)

The two-dimensional Laplacian in polar coordinates, ∇² = ∂²/∂r² + (1/r) ∂/∂r + (1/r²) ∂²/∂θ² (Sadd 7.6.5).
"""
laplacian_polar(f) = expand_derivatives(Dr(Dr(f)) + Dr(f)/r + Dθ(Dθ(f))/r^2)

# ╔═╡ d4000000-0000-0000-0000-000000000008
"""
    biharmonic_polar(ϕ)

∇⁴ϕ = ∇²(∇²ϕ) in polar coordinates (Sadd 7.6.8).
"""
biharmonic_polar(ϕ) = simplify(laplacian_polar(laplacian_polar(ϕ)))

# ╔═╡ d4000000-0000-0000-0000-000000000009
"""
    stresses_polar(ϕ) -> (σr, σθ, τrθ)

Stresses from the Airy stress function in polar coordinates (Sadd 7.6.7).
"""
stresses_polar(ϕ) = (simplify(expand_derivatives(Dr(ϕ)/r + Dθ(Dθ(ϕ))/r^2)),
                     simplify(expand_derivatives(Dr(Dr(ϕ)))),
                     simplify(expand_derivatives(-Dr(Dθ(ϕ)/r))))

# ╔═╡ d4000000-0000-0000-0000-00000000000a
md"""
Helpers. `coeffs_in` is the same as in Lecture 12: a condition that must hold *for all* ``r`` is equivalent to every coefficient of its polynomial in ``r`` vanishing. The integration helper changes: the stresses in polar coordinates are not polynomials in ``r`` but **Laurent polynomials** (powers ``r^{-2}, r^0, r^2, \dots``), possibly multiplied by ``\ln r``. `∫r` clears the negative powers by multiplying through by ``r^k``, reads off the coefficients, and integrates term by term using ``\int r^m dr = r^{m+1}/(m+1)``, ``\int r^{-1} dr = \ln r`` and ``\int r^m \ln r\, dr = \frac{r^{m+1}}{m+1}\ln r - \frac{r^{m+1}}{(m+1)^2}``.
"""

# ╔═╡ d4000000-0000-0000-0000-00000000000b
begin
	coeffs_in(expr, vars...) = collect(values(first(Symbolics.polynomial_coeffs(expand(expr), collect(vars)))))

	function ∫r(expr; k = 4)
		ℓ = only(@variables ℓ)                                                # stands in for log(r) while we read coefficients
		e = expand(Symbolics.simplify_fractions(expand(expr) * r^k))
		e = expand(substitute(e, Dict(log(r) => ℓ)))
		d, rem = Symbolics.polynomial_coeffs(e, [r, ℓ])
		isequal(Symbolics.value(rem), 0) || error("expression is not a Laurent polynomial in r (times log r)")
		s = 0
		for (mono, co) in d
			m = Symbolics.degree(mono, r) - k
			n = Symbolics.degree(mono, ℓ)
			n ≤ 1 || error("powers of log r above the first are not handled")
			if n == 0
				s += m == -1 ? co*log(r) : co*r^(m + 1)/(m + 1)
			else
				s += m == -1 ? co*log(r)^2/2 : co*(r^(m + 1)/(m + 1)*log(r) - r^(m + 1)/(m + 1)^2)
			end
		end
		return expand(s)
	end

	∫r(expr, a, b; k = 4) = expand(substitute(∫r(expr; k), Dict(r => b)) - substitute(∫r(expr; k), Dict(r => a)))
end

# ╔═╡ d4000000-0000-0000-0000-00000000000c
"""
    solve_linear_bcs(eqs, unknowns) -> (solution::Dict, n_redundant)

Exact Gauss–Jordan elimination on the augmented matrix [A | b] of a linear system that may have
more equations than unknowns (same as Lecture 12). Numeric pivots are preferred, so we never divide by a
symbolic parameter unless we must. Errors if the system is inconsistent or leaves an unknown undetermined.
"""
function solve_linear_bcs(eqs, unknowns)
	isz(a)  = isequal(Symbolics.value(simplify(a)), 0)
	isnum(a) = Symbolics.value(a) isa Number
	lhs = [eq.lhs - eq.rhs for eq in eqs]             # every entry must equal zero
	A   = Symbolics.jacobian(lhs, unknowns)            # linear in the unknowns ⇒ A is free of them
	b   = -substitute.(lhs, Ref(Dict(unknowns .=> 0)))
	M   = Num.(hcat(A, b))
	m, n = size(A)
	row, pivots = 1, Int[]
	for col in 1:n
		cand = [i for i in row:m if !isz(M[i, col])]
		isempty(cand) && continue
		k = findfirst(i -> isnum(M[i, col]), cand)
		i₀ = k === nothing ? cand[1] : cand[k]
		M[row, :], M[i₀, :] = M[i₀, :], M[row, :]
		M[row, :] = simplify.(M[row, :] ./ M[row, col])
		for i in 1:m
			i == row && continue
			M[i, :] = simplify.(M[i, :] .- M[i, col] .* M[row, :])
		end
		push!(pivots, col)
		row += 1
	end
	all(isz, M[row:end, end]) || error("the boundary conditions are inconsistent")
	length(pivots) == n || error("$(n - length(pivots)) coefficient(s) left undetermined")
	return Dict(unknowns[pivots[i]] => M[i, end] for i in eachindex(pivots)), m - length(pivots)
end

# ╔═╡ d4000000-0000-0000-0000-00000000000d
md"""
## 4. Trial stress function: the axisymmetric terms of the Michell solution

The Michell solution (Sadd 8.3.6) is the complete list of separable solutions of the polar biharmonic equation. The terms whose **stresses** do not depend on ``\theta`` are

```math
\phi = a_0 + a_1\ln r + a_2 r^2 + a_3 r^2 \ln r + a_4\,\theta
```

``a_0`` produces no stress and is dropped. The last term deserves a comment: ``\phi = a_4\theta`` is not single-valued on a full ring, but its stresses are (``\tau_{r\theta} = a_4/r^2``, a pure axisymmetric shear — Sadd §8.3.2), and on our half ring ``\theta`` is a perfectly good function. Sadd solves Example 8-6 having already discarded ``a_3`` and ``a_4``: ``a_4`` because there is no shear load, ``a_3`` because on a *full* ring its tangential displacement is multivalued (it grows by ``4\pi a_3 r/E`` every time round the ring). Our domain does **not** encircle the origin, so neither argument is available to us a priori. We keep all four terms and let the boundary conditions decide.
"""

# ╔═╡ d4000000-0000-0000-0000-00000000000e
begin
	@variables a1 a2 a3 a4
	ϕ_trial = a1*log(r) + a2*r^2 + a3*r^2*log(r) + a4*θ
end

# ╔═╡ d4000000-0000-0000-0000-00000000000f
md"**Biharmonic check.** Every term of the Michell solution is biharmonic by construction, so this should come out zero without any condition on the coefficients:"

# ╔═╡ d4000000-0000-0000-0000-000000000010
∇⁴ϕ = biharmonic_polar(ϕ_trial)

# ╔═╡ d4000000-0000-0000-0000-000000000011
md"**Stresses from the trial function** — Sadd (8.3.8) plus the ``a_4`` shear:"

# ╔═╡ d4000000-0000-0000-0000-000000000012
σr_trial, σθ_trial, τrθ_trial = stresses_polar(ϕ_trial)

# ╔═╡ d4000000-0000-0000-0000-000000000013
md"""
```math
\sigma_r = \frac{a_1}{r^2} + 2a_2 + a_3\left(1 + 2\ln r\right), \qquad
\sigma_\theta = -\frac{a_1}{r^2} + 2a_2 + a_3\left(3 + 2\ln r\right), \qquad
\tau_{r\theta} = \frac{a_4}{r^2}
```

**Plane strain strains** of the trial field:
"""

# ╔═╡ d4000000-0000-0000-0000-000000000014
begin
	εr_trial  = simplify((1 + ν)/E * ((1 - ν)*σr_trial - ν*σθ_trial))
	εθ_trial  = simplify((1 + ν)/E * ((1 - ν)*σθ_trial - ν*σr_trial))
	γrθ_trial = simplify(2(1 + ν)/E * τrθ_trial)
	(εr_trial, εθ_trial, γrθ_trial)
end

# ╔═╡ d4000000-0000-0000-0000-000000000015
md"""
## 5. Displacements of the general axisymmetric field

Before imposing any boundary condition, integrate the strain–displacement relations for the trial field. The algebra is short enough to follow by hand, with `Symbolics.jl` doing the bookkeeping.

**Step 1 — radial strain.** ``\partial u_r/\partial r = \varepsilon_r(r)`` gives ``u_r = \int \varepsilon_r\,dr + f(\theta)``, where the "constant" of integration may depend on ``\theta``:
"""

# ╔═╡ d4000000-0000-0000-0000-000000000016
ur_part = ∫r(εr_trial)

# ╔═╡ d4000000-0000-0000-0000-000000000017
md"""
**Step 2 — hoop strain.** ``u_r + \partial u_\theta/\partial\theta = r\,\varepsilon_\theta``, so ``\partial u_\theta/\partial \theta = r\varepsilon_\theta - \int\varepsilon_r\,dr - f(\theta)``. The combination ``r\varepsilon_\theta - \int\varepsilon_r\,dr`` turns out to be remarkably simple — everything except the ``a_3`` term cancels:
"""

# ╔═╡ d4000000-0000-0000-0000-000000000018
κr = simplify(expand(r*εθ_trial - ur_part))

# ╔═╡ d4000000-0000-0000-0000-000000000019
md"""
so ``\partial u_\theta/\partial\theta = \kappa r - f(\theta)`` with ``\kappa = 4(1-\nu^2)a_3/E``, and ``u_\theta = \kappa\, r\,\theta - F(\theta) + g(r)`` with ``F' = f``.

**Step 3 — shear strain.** Substitute both displacements into ``\gamma_{r\theta} = \frac{1}{r}\partial_\theta u_r + \partial_r u_\theta - u_\theta/r``. The ``\kappa r\theta`` terms cancel (``\partial_r(\kappa r\theta) = \kappa r\theta / r``) and what is left is

```math
\frac{f'(\theta) + F(\theta)}{r} + g'(r) - \frac{g(r)}{r} = \gamma_{r\theta} = \frac{2(1+\nu)\,a_4}{E\,r^2}
```

The left side separates into a function of ``\theta`` and a function of ``r``, so each must be constant. ``F'' + F = \text{const}`` gives ``f = u_{x0}\cos\theta + u_{y0}\sin\theta``: a rigid translation ``(u_{x0}, u_{y0})`` written in polar components. ``g' - g/r = 2(1+\nu)a_4/(Er^2) - \text{const}/r`` gives ``g = -(1+\nu)a_4/(Er) + \omega_0 r``: the displacement of the ``a_4`` shear field plus a rigid rotation ``\omega_0``. Altogether,

```math
u_r = \int\varepsilon_r\,dr + u_{x0}\cos\theta + u_{y0}\sin\theta, \qquad
u_\theta = \frac{4(1-\nu^2)a_3}{E}\,r\theta - \frac{(1+\nu)a_4}{E\,r} - u_{x0}\sin\theta + u_{y0}\cos\theta + \omega_0 r
```

This is Sadd (8.3.9) converted to plane strain with Table 7-1, plus the ``a_4`` term. Assemble it and let `Symbolics` confirm that it reproduces all three strains (each residual should be zero):
"""

# ╔═╡ d4000000-0000-0000-0000-00000000001a
begin
	@variables ux0 uy0 ω0
	ur_general = ur_part + ux0*cos(θ) + uy0*sin(θ)
	uθ_general = 4(1 - ν^2)*a3*r*θ/E - (1 + ν)*a4/(E*r) - ux0*sin(θ) + uy0*cos(θ) + ω0*r
	(ur_general, uθ_general)
end

# ╔═╡ d4000000-0000-0000-0000-00000000001b
strain_check_general = simplify.(expand.([
	expand_derivatives(Dr(ur_general)) - εr_trial,
	expand_derivatives((ur_general + Dθ(uθ_general))/r) - εθ_trial,
	expand_derivatives(Dθ(ur_general)/r + Dr(uθ_general) - uθ_general/r) - γrθ_trial,
]))

# ╔═╡ d4000000-0000-0000-0000-00000000001c
md"""
## 6. What the supports say

Impose the support and symmetry conditions on the general field. Each must hold **for every** ``r`` across the foot, so we multiply through by ``r`` (or ``r^2``) to clear the ``1/r`` terms and match coefficients of powers of ``r``:

* right foot: ``u_\theta(r, 0) = 0`` and ``\tau_{r\theta}(r, 0) = 0``;
* left foot: ``u_\theta(r, \pi) = 0`` and ``\tau_{r\theta}(r, \pi) = 0``;
* axis of symmetry: ``u_\theta(r, \pi/2) = 0``.

(`at_angle` substitutes the exact values of ``\sin\theta`` and ``\cos\theta`` before replacing ``\theta`` itself, so no floating-point ``\sin\pi = 10^{-16}`` noise enters the exact row reduction. ``\pi`` is kept as a symbol.)
"""

# ╔═╡ d4000000-0000-0000-0000-00000000001d
begin
	at_angle(expr, θ₀, sinθ₀, cosθ₀) = substitute(substitute(expr, Dict(sin(θ) => sinθ₀, cos(θ) => cosθ₀)), Dict(θ => θ₀))
	π_ = Num(π)

	support_eqs = vcat(
		# right foot θ = 0, for all r:  uθ = 0,  τrθ = 0
		[e ~ 0 for e in coeffs_in(r   * at_angle(uθ_general, 0, 0, 1), r)],
		[e ~ 0 for e in coeffs_in(r^2 * at_angle(τrθ_trial,  0, 0, 1), r)],
		# left foot θ = π, for all r:  uθ = 0,  τrθ = 0
		[e ~ 0 for e in coeffs_in(r   * at_angle(uθ_general, π_, 0, -1), r)],
		[e ~ 0 for e in coeffs_in(r^2 * at_angle(τrθ_trial,  π_, 0, -1), r)],
		# axis of symmetry θ = π/2, for all r:  uθ = 0
		[e ~ 0 for e in coeffs_in(r   * at_angle(uθ_general, π_/2, 1, 0), r)],
	)
end

# ╔═╡ d4000000-0000-0000-0000-00000000001e
support_map, n_redundant_support = solve_linear_bcs(support_eqs, [a3, a4, ux0, uy0, ω0])

# ╔═╡ d4000000-0000-0000-0000-00000000001f
md"""
$(length(support_eqs)) equations, $(n_redundant_support) of them redundant, and every one of the five unknowns is zero. The rigid-body motion is gone, ``a_4 = 0`` because the feet carry no shear, and — the interesting one — ``a_3 = 0``.

The ``r^2\ln r`` term is killed not by single-valuedness, as it is in the full ring, but by the **supports**: its tangential displacement ``4(1-\nu^2)a_3 r\theta/E`` grows linearly with ``\theta``, so if the right foot is held at its level the left foot would rise or drop by ``4\pi(1-\nu^2)a_3 r/E``. Physically, ``a_3 \neq 0`` is the self-stress of an arch that has been forced open or closed by a prescribed angle at the feet (Sadd Exercise 8-21, the "dislocation" solution) — a stress field that exists with no pressure at all. Both feet held at the same level means no such opening, hence ``a_3 = 0``. Had the supports been at different levels, or had one foot been forced to rotate, this term would be exactly what we needed.
"""

# ╔═╡ d4000000-0000-0000-0000-000000000020
md"""
## 7. What the pressure says

With ``a_3 = a_4 = 0`` only the Lamé terms ``a_1\ln r + a_2 r^2`` survive. Two pressure conditions for two unknowns; the shear conditions on both surfaces are satisfied identically:
"""

# ╔═╡ d4000000-0000-0000-0000-000000000021
begin
	σr_1  = substitute(σr_trial,  support_map)
	τrθ_1 = substitute(τrθ_trial, support_map)
	pressure_eqs = [
		substitute(σr_1,  Dict(r => R))     ~ 0,      # inner surface free
		substitute(σr_1,  Dict(r => R + t)) ~ -p,     # ground pressure on the outer surface
		substitute(τrθ_1, Dict(r => R))     ~ 0,
		substitute(τrθ_1, Dict(r => R + t)) ~ 0,
	]
end

# ╔═╡ d4000000-0000-0000-0000-000000000022
pressure_map, n_redundant_pressure = solve_linear_bcs(pressure_eqs, [a1, a2])

# ╔═╡ d4000000-0000-0000-0000-000000000023
md"""
```math
a_1 = \frac{p\,R^2(R+t)^2}{t\,(2R + t)}, \qquad a_2 = -\frac{p\,(R+t)^2}{2\,t\,(2R + t)}
```

(``t(2R+t) = (R+t)^2 - R^2``; these are Sadd's ``A = a_1`` and ``B = 2a_2`` of (8.4.2) with ``r_1 = R``, ``r_2 = R + t``, ``p_1 = 0``, ``p_2 = p``.) The stress function and the stress field:
"""

# ╔═╡ d4000000-0000-0000-0000-000000000024
begin
	coeff_map = merge(support_map, pressure_map)
	ϕ = simplify(substitute(ϕ_trial, coeff_map))
end

# ╔═╡ d4000000-0000-0000-0000-000000000025
begin
	σr, σθ, τrθ = stresses_polar(ϕ)
	σz = simplify(ν*(σr + σθ))             # plane strain: the stress that keeps εz = 0
	(σr = σr, σθ = σθ, τrθ = τrθ, σz = σz)
end

# ╔═╡ d4000000-0000-0000-0000-000000000026
md"""
Written with the common factor ``K = \dfrac{p\,(R+t)^2}{t\,(2R+t)}``:

```math
\sigma_r = -K\left(1 - \frac{R^2}{r^2}\right), \qquad
\sigma_\theta = -K\left(1 + \frac{R^2}{r^2}\right), \qquad
\tau_{r\theta} = 0, \qquad
\sigma_z = -2\nu K
```

Sadd (8.4.3)–(8.4.4) with the roles of the two pressures reversed from his Figure 8-9. Observations:

* **Both in-plane stresses are compressive everywhere**, and ``\sigma_\theta`` dominates: the pressure ``p`` is small compared with the hoop stress it produces, by the factor ``\sim R/t`` familiar from thin rings.
* The **largest** stress is the hoop stress at the **inner** surface, ``\sigma_\theta(R) = -2K``, on the surface that carries no load at all. The outer surface, where the pressure acts, has the smaller ``|\sigma_\theta(R+t)| = K(1 + R^2/(R+t)^2)``.
* ``\sigma_z`` is uniform through the thickness and is the only stress that knows about the material (through ``\nu``). It is compressive: the tunnel would like to lengthen under the Poisson effect of the in-plane compression, and ``\sigma_z`` stops it.
* The stresses split into two classic fields: the ``-K`` part is uniform hydrostatic compression (the ``a_2 r^2`` term, Sadd §8.4.2 in reverse), and the ``\pm K R^2/r^2`` part is the *pressurised hole* field of Sadd §8.4.1 (the ``a_1 \ln r`` term) — pure shear in the ``r``–``\theta`` plane, which is why its contributions to ``\sigma_r`` and ``\sigma_\theta`` are equal and opposite.

Verify every condition we imposed (each residual should be zero):
"""

# ╔═╡ d4000000-0000-0000-0000-000000000027
bc_residuals = simplify.(expand.([
	biharmonic_polar(ϕ),
	substitute(σr,  Dict(r => R)),
	substitute(σr,  Dict(r => R + t)) + p,
	τrθ,                                                   # for all r, θ — covers both surfaces and both feet
	substitute(ν*(σr + σθ) - σz, Dict(r => R)),            # σz is the plane strain value
]))

# ╔═╡ d4000000-0000-0000-0000-000000000028
md"""
**Do the supports carry the load?** We never prescribed the reaction. The vertical force per unit length of tunnel that each foot delivers is the integral of the hoop stress across the foot, ``N = \int_R^{R+t}\sigma_\theta(r, 0)\,dr``. Statics of the half ring says the total downward load is the vertical resultant of a radial pressure ``p`` acting on a semicircle of radius ``R + t``, which is ``2p(R+t)``, so each foot must carry ``p(R+t)``:
"""

# ╔═╡ d4000000-0000-0000-0000-000000000029
begin
	N_foot = simplify(∫r(σθ, R, R + t))                               # hoop force per unit length (compression negative)
	M_foot = simplify(∫r(σθ*(r - (R + t/2)), R, R + t))              # its moment about mid-thickness
	(N_foot = N_foot, equilibrium_residual = simplify(N_foot + p*(R + t)), M_foot = M_foot)
end

# ╔═╡ d4000000-0000-0000-0000-00000000002a
md"""
``N = -p(R+t)`` exactly: equilibrium is built into the Airy stress function, so the reaction took care of itself. The moment of the hoop stress about the mid-thickness is not zero, because ``\sigma_\theta`` is larger at the inner face — the resultant acts slightly *inside* mid-thickness. The thrust in the arch is therefore eccentric, by a small amount that a thin-ring analysis cannot see.

## 8. Strains and displacements

Plane strain Hooke's law for the final stresses, then the displacements from the general field of Section 5 with every coefficient substituted:
"""

# ╔═╡ d4000000-0000-0000-0000-00000000002b
begin
	εr  = simplify((1 + ν)/E * ((1 - ν)*σr - ν*σθ))
	εθ  = simplify((1 + ν)/E * ((1 - ν)*σθ - ν*σr))
	γrθ = simplify(2(1 + ν)/E * τrθ)
	εz  = 0                                   # by definition of plane strain
	(εr = εr, εθ = εθ, γrθ = γrθ, εz = εz)
end

# ╔═╡ d4000000-0000-0000-0000-00000000002c
begin
	ur = simplify(substitute(ur_general, coeff_map))
	uθ = simplify(substitute(uθ_general, coeff_map))
	(ur = ur, uθ = uθ)
end

# ╔═╡ d4000000-0000-0000-0000-00000000002d
md"""
```math
u_r = -\frac{(1+\nu)\,K}{E}\left[(1-2\nu)\,r + \frac{R^2}{r}\right], \qquad u_\theta = 0, \qquad u_z = 0
```

(Sadd 8.4.5 with ``A = a_1``, ``B = 2a_2``.) Every point moves **radially toward the centre**: the crown comes down, the haunches move in, and the feet slide horizontally inward on their smooth supports by ``|u_r(r, 0)|`` — the sliding our support idealisation allows. Tunnel engineers call the resulting decrease in diameter, ``2|u_r(R)|``, the *convergence* of the lining.

The two terms are the displacements of the two stress fields identified above: ``(1-2\nu)r`` is the uniform contraction under hydrostatic compression, with the plane strain bulk factor ``(1+\nu)(1-2\nu)/E`` (it vanishes for an incompressible material, ``\nu = 1/2``, as it must); ``R^2/r`` is the volume-preserving ``1/r`` displacement of the pressurised-hole field (Sadd 8.4.8).

Check that the displacements reproduce the strains:
"""

# ╔═╡ d4000000-0000-0000-0000-00000000002e
strain_check = simplify.(expand.([
	expand_derivatives(Dr(ur)) - εr,
	expand_derivatives((ur + Dθ(uθ))/r) - εθ,
	expand_derivatives(Dθ(ur)/r + Dr(uθ) - uθ/r) - γrθ,
]))

# ╔═╡ d4000000-0000-0000-0000-00000000002f
md"""
**Plane stress, for comparison.** Had we (wrongly) treated the tunnel as a thin slice in plane stress, the stresses ``\sigma_r, \sigma_\theta`` would be identical, ``\sigma_z`` would be zero, and the radial displacement would be Sadd (8.3.9) with ``a_3 = 0``:

```math
u_r^{\text{plane stress}} = \frac{1}{E}\left[2(1-\nu)\,a_2\,r - (1+\nu)\frac{a_1}{r}\right]
```

The plane strain constraint ``\varepsilon_z = 0`` makes the lining stiffer in its own plane, so ``|u_r|`` is smaller in plane strain. We will quantify the difference in the numerical example.
"""

# ╔═╡ d4000000-0000-0000-0000-000000000030
ur_plane_stress = simplify(substitute((2(1 - ν)*a2*r - (1 + ν)*a1/r)/E, coeff_map))

# ╔═╡ d4000000-0000-0000-0000-000000000031
md"""
## 9. Comparison with the thin-ring formula

Strength of materials treats the lining as a thin ring: the uniform external pressure ``p`` acting at radius ``R + t`` produces a hoop compression ``N = -p(R+t)`` per unit length (the same equilibrium we just verified), spread *uniformly* over the thickness,

```math
\sigma_\theta^{\text{thin}} = \frac{N}{t} = -\frac{p\,(R+t)}{t}
```

and, neglecting ``\sigma_r`` (it is at most ``p``, small compared with ``\sigma_\theta``), the plane strain hoop strain ``\varepsilon_\theta = (1-\nu^2)\sigma_\theta/E`` at the mean radius ``r_m = R + t/2`` gives a radial displacement

```math
u_r^{\text{thin}} = r_m\,\varepsilon_\theta = -\frac{(1-\nu^2)\,p\,(R+t)\,(R + t/2)}{E\,t}
```

Form the differences symbolically. The *average* of the exact hoop stress through the thickness is ``N/t`` exactly — the thin-ring stress is the mean of the true one — and the exact stress at the inner and outer faces straddles it:
"""

# ╔═╡ d4000000-0000-0000-0000-000000000032
begin
	σθ_thin = -p*(R + t)/t
	ur_thin = -(1 - ν^2)*p*(R + t)*(R + t/2)/(E*t)
	thin_ring = (
		mean_exact_minus_thin = simplify(N_foot/t - σθ_thin),
		inner_minus_thin      = simplify(substitute(σθ, Dict(r => R)) - σθ_thin),
		outer_minus_thin      = simplify(substitute(σθ, Dict(r => R + t)) - σθ_thin),
		ratio_inner_to_outer  = simplify(substitute(σθ, Dict(r => R)) / substitute(σθ, Dict(r => R + t))),
	)
end

# ╔═╡ d4000000-0000-0000-0000-00000000003e
md"""
```math
\sigma_\theta(R) - \sigma_\theta^{\text{thin}} = -\frac{p\,(R+t)}{2R+t}, \qquad
\sigma_\theta(R+t) - \sigma_\theta^{\text{thin}} = +\frac{p\,R}{2R+t}
```

Both deviations are of the order of the pressure ``p`` itself, i.e. smaller than the hoop stress by a factor ``\sim t/(2R)`` — the thin-ring formula is good to a few percent whenever the lining is thin. The inner-to-outer ratio depends only on ``t/R``:

```math
\frac{\sigma_\theta(R)}{\sigma_\theta(R+t)} = \frac{2(R+t)^2}{(R+t)^2 + R^2}
```

For a thin lining it tends to 1 — uniform hoop stress, the thin-ring picture. For a very thick one (``t \gg R``) it tends to 2: the stress concentration factor of a hole in an infinite medium under all-round pressure (Sadd 8.4.10 with the sign of ``T`` reversed). The lining in the example below sits much closer to the first limit.
"""

# ╔═╡ d4000000-0000-0000-0000-000000000033
md"""
## 10. Numerical example

A cast-in-place **concrete** lining: inner radius ``R = 3`` m, thickness ``t = 300`` mm, ``E = 30`` GPa, ``\nu = 0.2``. The ground pressure is taken as ``p = 0.5`` MPa — roughly 25 m of soil at 20 kN/m³, applied uniformly all round. Units: N, mm, MPa; all forces are per mm of tunnel length.
"""

# ╔═╡ d4000000-0000-0000-0000-000000000034
prm = (p = 0.5,          # MPa  ground pressure
       R = 3000.0,       # mm   inner radius
       t = 300.0,        # mm   lining thickness
       E = 30_000.0,     # MPa  concrete
       ν = 0.2)

# ╔═╡ d4000000-0000-0000-0000-000000000035
begin
	# build_function compiles each symbolic expression to a callable f(r, p, R, t, E, ν)
	args = (r, p, R, t, E, ν)
	σr_f   = build_function(σr,   args...; expression = Val(false))
	σθ_f   = build_function(σθ,   args...; expression = Val(false))
	σz_f   = build_function(σz,   args...; expression = Val(false))
	εr_f   = build_function(εr,   args...; expression = Val(false))
	εθ_f   = build_function(εθ,   args...; expression = Val(false))
	ur_f   = build_function(ur,   args...; expression = Val(false))
	urps_f = build_function(ur_plane_stress, args...; expression = Val(false))
	ev(f, ri) = f(ri, prm.p, prm.R, prm.t, prm.E, prm.ν)
	Sr(ri) = ev(σr_f, ri);  Sθ(ri) = ev(σθ_f, ri);  Sz(ri) = ev(σz_f, ri)
	Er(ri) = ev(εr_f, ri);  Eθ(ri) = ev(εθ_f, ri)
	Ur(ri) = ev(ur_f, ri);  Urps(ri) = ev(urps_f, ri)
	σθ_thin_num = -prm.p*(prm.R + prm.t)/prm.t
	ur_thin_num = -(1 - prm.ν^2)*prm.p*(prm.R + prm.t)*(prm.R + prm.t/2)/(prm.E*prm.t)
	K_num = prm.p*(prm.R + prm.t)^2/(prm.t*(2prm.R + prm.t))
	(K = K_num, σθ_thin = σθ_thin_num, ur_thin = ur_thin_num)
end

# ╔═╡ d4000000-0000-0000-0000-000000000036
md"""
### Exact solution vs thin ring, side by side

``r_i = R`` is the inner (free) surface, ``r_o = R + t`` the outer (loaded) one, ``r_m`` the mid-thickness.
"""

# ╔═╡ d4000000-0000-0000-0000-000000000037
comparison = let
	ri, rm, ro = prm.R, prm.R + prm.t/2, prm.R + prm.t
	rows = [
		("σθ inner surface rᵢ  [MPa]",            Sθ(ri),            σθ_thin_num),
		("σθ mid-thickness rₘ  [MPa]",            Sθ(rm),            σθ_thin_num),
		("σθ outer surface rₒ  [MPa]",            Sθ(ro),            σθ_thin_num),
		("σθ average through thickness  [MPa]",   -prm.p*ro/prm.t,   σθ_thin_num),
		("σr inner surface rᵢ  [MPa]",            Sr(ri),            0.0),
		("σr outer surface rₒ  [MPa]",            Sr(ro),            0.0),
		("σz (uniform)  [MPa]",                   Sz(ri),            prm.ν*σθ_thin_num),
		("uᵣ inner surface rᵢ  [mm]",             Ur(ri),            ur_thin_num),
		("uᵣ outer surface rₒ  [mm]",             Ur(ro),            ur_thin_num),
		("uᵣ inner surface, plane stress  [mm]",  Urps(ri),          NaN),
		("convergence (decrease in inner diameter)  [mm]", 2abs(Ur(ri)), 2abs(ur_thin_num)),
		("εθ inner surface  [10⁻⁶]",              1e6*Eθ(ri),        1e6*(1 - prm.ν^2)*σθ_thin_num/prm.E),
		("εr inner surface  [10⁻⁶]",              1e6*Er(ri),        1e6*(-prm.ν*(1 + prm.ν))*σθ_thin_num/prm.E),
	]
	fmt(v) = isnan(v) ? "—" : string(round(v, digits = 4) + 0.0)     # + 0.0 turns -0.0 into 0.0
	hdr = "| quantity | elasticity (plane strain) | thin ring |\n|---|---|---|\n"
	Markdown.parse(hdr * join(["| $(q) | $(fmt(a)) | $(fmt(b)) |" for (q, a, b) in rows], "\n"))
end

# ╔═╡ d4000000-0000-0000-0000-000000000038
md"""
With ``t/R = `` $(prm.t/prm.R) the thin-ring hoop stress is the exact average, but it underestimates the inner-face stress by $(round(100*(abs(Sθ(prm.R)/σθ_thin_num) - 1), digits = 1)) % and overestimates the outer-face stress by $(round(100*(1 - abs(Sθ(prm.R + prm.t)/σθ_thin_num)), digits = 1)) %. The radial stress never exceeds ``p`` = $(prm.p) MPa, about a tenth of the hoop stress, which is why the thin-ring theory can ignore it. The displacement predicted by the thin ring is within $(round(100*abs(ur_thin_num/Ur(prm.R) - 1), digits = 1)) % of the exact inner-surface value. Plane stress would overestimate the convergence by $(round(100*(Urps(prm.R)/Ur(prm.R) - 1), digits = 1)) % — the ``\varepsilon_z = 0`` constraint matters even for ``\nu = 0.2``. The concrete stress of about $(round(abs(Sθ(prm.R)), digits = 1)) MPa is comfortably below a typical compressive strength, and the convergence of $(round(2abs(Ur(prm.R)), digits = 2)) mm on a 6 m tunnel is imperceptible.
"""

# ╔═╡ d4000000-0000-0000-0000-000000000039
md"### Stresses and displacement through the thickness"

# ╔═╡ d4000000-0000-0000-0000-00000000003a
fig_through = let
	rs = range(prm.R, prm.R + prm.t, length = 61)
	fig = Figure(size = (1000, 400))

	ax1 = Axis(fig[1, 1], title = "Stresses through the lining thickness",
	           xlabel = "r [mm]", ylabel = "stress [MPa]")
	lines!(ax1, rs, Sθ.(rs), color = :firebrick, linewidth = 3, label = "σθ")
	lines!(ax1, rs, Sr.(rs), color = :steelblue, linewidth = 3, label = "σr")
	lines!(ax1, rs, Sz.(rs), color = :darkgreen, linewidth = 3, label = "σz (plane strain)")
	hlines!(ax1, [σθ_thin_num], color = :black, linestyle = :dash, linewidth = 2, label = "σθ thin ring  −p(R+t)/t")
	hlines!(ax1, [0.0], color = (:gray, 0.5))
	axislegend(ax1, position = :rc)

	ax2 = Axis(fig[1, 2], title = "Radial displacement through the thickness",
	           xlabel = "r [mm]", ylabel = "uᵣ [mm]")
	lines!(ax2, rs, Ur.(rs), color = :firebrick, linewidth = 3, label = "uᵣ plane strain")
	lines!(ax2, rs, Urps.(rs), color = :firebrick, linewidth = 2, linestyle = :dot, label = "uᵣ plane stress")
	hlines!(ax2, [ur_thin_num], color = :black, linestyle = :dash, linewidth = 2, label = "uᵣ thin ring")
	axislegend(ax2, position = :rb)
	fig
end

# ╔═╡ d4000000-0000-0000-0000-00000000003b
md"""
The hoop stress is largest in magnitude on the free inner face and decays as ``1/r^2`` toward the outer face; the radial stress grows from zero to ``-p``; ``\sigma_z`` is flat. The displacement is nearly uniform across the thickness (``|u_r|`` is slightly larger inside because the inner fibre is more compressed hoop-wise), and plane strain is visibly stiffer than plane stress.

### Hoop stress over the cross-section and the deformed shape

Left: ``\sigma_\theta`` painted on the arch — axisymmetric, so it varies only with ``r``. Right: the deformed lining with displacements exaggerated, together with the mirror image that completes the full ring of Section 2. Every point moves radially inward; the feet slide horizontally on their supports and stay at their level.
"""

# ╔═╡ d4000000-0000-0000-0000-00000000003c
fig_section = let
	fig = Figure(size = (1100, 480))

	# --- hoop stress painted on the arch, as a patchwork of polar quads ---
	ax1 = Axis(fig[1, 1], title = "σθ [MPa] over the cross-section", xlabel = "x [mm]", ylabel = "y [mm]", aspect = DataAspect())
	nr, nθ = 10, 90
	redges = range(prm.R, prm.R + prm.t, length = nr + 1)
	θedges = range(0, π, length = nθ + 1)
	quads, vals = Vector{Point2f}[], Float64[]
	for i in 1:nr, j in 1:nθ
		ra, rb = redges[i], redges[i + 1]
		θa, θb = θedges[j], θedges[j + 1]
		push!(quads, Point2f[(ra*cos(θa), ra*sin(θa)), (rb*cos(θa), rb*sin(θa)), (rb*cos(θb), rb*sin(θb)), (ra*cos(θb), ra*sin(θb))])
		push!(vals, Sθ((ra + rb)/2))
	end
	poly!(ax1, quads, color = vals, colormap = :viridis, strokewidth = 0)
	Colorbar(fig[1, 2], limits = extrema(vals), colormap = :viridis, label = "σθ [MPa]")
	# pressure arrows on the outer surface
	for θa in range(0.05π, 0.95π, length = 13)
		ro = prm.R + prm.t
		arrows2d!(ax1, [1.2ro*cos(θa)], [1.2ro*sin(θa)], [-0.18ro*cos(θa)], [-0.18ro*sin(θa)], color = :seagreen)
	end
	xlims!(ax1, -1.3(prm.R + prm.t), 1.3(prm.R + prm.t)); ylims!(ax1, -0.1prm.R, 1.3(prm.R + prm.t))

	# --- deformed shape, exaggerated ---
	scale = 1500
	ax2 = Axis(fig[1, 3], title = "Deformed shape, displacements × $scale", xlabel = "x [mm]", ylabel = "y [mm]", aspect = DataAspect())
	θs = range(0, π, length = 181)
	arc(rad, s, sgn) = [Point2f((rad + s*Ur(rad))*cos(θi), sgn*(rad + s*Ur(rad))*sin(θi)) for θi in θs]
	foot(s, sgn) = [Point2f(sgn*(prm.R + s*Ur(prm.R)), 0), Point2f(sgn*(prm.R + prm.t + s*Ur(prm.R + prm.t)), 0)]
	for (s, col, lw, lab) in ((0, (:gray, 0.8), 2, "undeformed"), (scale, :firebrick, 3, "deformed"))
		lines!(ax2, arc(prm.R, s, 1), color = col, linewidth = lw, label = lab)
		lines!(ax2, arc(prm.R + prm.t, s, 1), color = col, linewidth = lw)
		lines!(ax2, foot(s, 1), color = col, linewidth = lw)
		lines!(ax2, foot(s, -1), color = col, linewidth = lw)
		# mirror image: the lower half of the full ring that the smooth supports stand in for
		lines!(ax2, arc(prm.R, s, -1), color = col, linewidth = 1, linestyle = :dash)
		lines!(ax2, arc(prm.R + prm.t, s, -1), color = col, linewidth = 1, linestyle = :dash)
	end
	hlines!(ax2, [0.0], color = :black, linewidth = 1)   # the support line
	axislegend(ax2, position = :rb)
	fig
end

# ╔═╡ d4000000-0000-0000-0000-00000000003f
md"""
## 11. Fixed, non-sliding supports

Suppose instead that the feet are **cast into a rigid invert**, or gripped by the rock, so that they can neither slide, nor lift, nor rotate. The support conditions of Section 2 become, pointwise across both feet,

```math
u_r(r, 0) = u_\theta(r, 0) = 0, \qquad u_r(r, \pi) = u_\theta(r, \pi) = 0, \qquad R \le r \le R + t
```

and the shear traction is no longer prescribed: the support supplies whatever ``\tau_{r\theta}`` it takes to hold the foot. Displacements are given on the feet and tractions on the curved surfaces — a *mixed* boundary value problem.

**Why this problem has no closed-form solution.** Every term of the Michell solution gives a foot displacement ``u_r(r, 0)`` that is a short Laurent polynomial in ``r`` (powers ``r^{-1}, r, r^3`` and ``r\ln r``). With finitely many terms it cannot vanish identically across the foot while the tractions also vanish on the two curved surfaces: the exact solution needs an infinite series of Michell terms, and at the two corners where a clamped face meets a free face it has the mild stress singularity of a clamped–free wedge (Williams 1952; Barber, *Elasticity*, §11.2) that no finite sum can reproduce.

**What we can do: fix the feet in the Saint-Venant sense.** Replace the pointwise conditions by their resultants. The foot has zero *mean* horizontal displacement, zero *mean* vertical displacement and zero *mean* rotation:

```math
\Delta = \frac{1}{t}\int_R^{R+t} u_r(r, 0)\,dr = 0, \qquad
\Lambda = \frac{1}{t}\int_R^{R+t} u_\theta(r, 0)\,dr = 0, \qquad
\Theta = \frac{12}{t^3}\int_R^{R+t} u_\theta(r, 0)\,(r - r_m)\,dr = 0, \qquad r_m = R + \tfrac{t}{2}
```

The factor ``12/t^3`` makes ``\Theta`` equal to the angle ``\omega`` when the foot rotates rigidly, ``u_\theta = \omega(r - r_m)``. This is exactly what Sadd does at the ends of every beam in Chapter 8 — Example 8-9 uses "resultant force boundary conditions at the ends and exact pointwise specifications along the lateral curved boundaries" — and it is what strength of materials means by a fixed end. By Saint-Venant's principle the solution obtained this way differs from the exact one only within about one thickness ``t`` of each foot. Everywhere else it *is* the answer.

Three conditions per foot instead of two means the axisymmetric terms will not be enough: the supports can now pull the feet horizontally and bend them, and the stresses will depend on ``\theta``.

### 11.1 Trial stress function: the symmetric ``\theta``-dependent Michell terms

The mirror symmetry about the crown ``\theta = \pi/2`` is intact, so we add the Michell terms (Sadd 8.3.6) with one angular wavelength whose *stresses* are symmetric about the crown — the ``\sin\theta`` family and the Flamant term ``r\theta\cos\theta``:

```math
\phi = a_1\ln r + a_2 r^2 + a_3 r^2\ln r + a_4\theta + \left(b_1 r^3 + \frac{b_2}{r} + b_3\, r\ln r\right)\sin\theta + b_4\, r\theta\cos\theta
```

``r\sin\theta = y`` produces no stress and is omitted; the ``\cos\theta`` family and ``r\theta\sin\theta`` are antisymmetric about the crown and would belong to a sideways load. We keep ``b_4`` to let the boundary conditions decide its fate. From here on ``\pi`` is carried as the symbol ``\Pi``, so that the exact row reduction never meets a floating-point number; it becomes ``3.14159\ldots`` only in the numerical example.
"""

# ╔═╡ d4000000-0000-0000-0000-000000000040
begin
	@variables b1 b2 b3 b4 Π
	ϕ_b = (b1*r^3 + b2/r + b3*r*log(r))*sin(θ) + b4*r*θ*cos(θ)
	ϕ_trial_fixed = ϕ_trial + ϕ_b
end

# ╔═╡ d4000000-0000-0000-0000-000000000041
∇⁴ϕ_fixed = biharmonic_polar(ϕ_trial_fixed)

# ╔═╡ d4000000-0000-0000-0000-000000000042
σr_tf, σθ_tf, τrθ_tf = stresses_polar(ϕ_trial_fixed)

# ╔═╡ d4000000-0000-0000-0000-000000000043
md"""
```math
\sigma_r = \frac{a_1}{r^2} + 2a_2 + a_3(1 + 2\ln r) + \left(2b_1 r - \frac{2b_2}{r^3} + \frac{b_3 - 2b_4}{r}\right)\sin\theta, \qquad
\sigma_\theta = -\frac{a_1}{r^2} + 2a_2 + a_3(3 + 2\ln r) + \left(6b_1 r + \frac{2b_2}{r^3} + \frac{b_3}{r}\right)\sin\theta
```
```math
\tau_{r\theta} = \frac{a_4}{r^2} - \left(2b_1 r - \frac{2b_2}{r^3} + \frac{b_3}{r}\right)\cos\theta
```

### 11.2 What the curved surfaces say

The pressure conditions hold **for every** ``\theta``, so on each surface the ``\theta``-independent part and the coefficient of ``\sin\theta`` (or ``\cos\theta``) must vanish separately. `trig` replaces ``\sin\theta, \cos\theta`` by placeholder symbols so that `coeffs_in` can read those coefficients off:
"""

# ╔═╡ d4000000-0000-0000-0000-000000000044
begin
	@variables s_ c_
	trig(expr) = substitute(expr, Dict(sin(θ) => s_, cos(θ) => c_))
	surface_eqs = vcat(
		[e ~ 0 for e in coeffs_in(trig(substitute(σr_tf,  Dict(r => R))), s_, c_)],          # inner surface free
		[e ~ 0 for e in coeffs_in(trig(substitute(σr_tf,  Dict(r => R + t))) + p, s_, c_)],  # ground pressure outside
		[e ~ 0 for e in coeffs_in(trig(substitute(τrθ_tf, Dict(r => R))), s_, c_)],
		[e ~ 0 for e in coeffs_in(trig(substitute(τrθ_tf, Dict(r => R + t))), s_, c_)],
	)
end

# ╔═╡ d4000000-0000-0000-0000-000000000045
surface_map, n_redundant_surface = solve_linear_bcs(surface_eqs, [a1, a2, a4, b1, b2, b4])

# ╔═╡ d4000000-0000-0000-0000-000000000046
md"""
$(length(surface_eqs)) equations, $(n_redundant_surface) redundant, six coefficients determined and **two left free**: ``a_3`` and ``b_3``. Writing ``a = R`` and ``b = R + t`` for the two radii:

* ``a_4 = b_4 = 0``: the axisymmetric shear ``a_4/r^2`` and Flamant's ``\sin\theta/r`` radial stress both put traction on the curved surfaces, so they are gone before the supports are even mentioned. (In Section 6 we let the supports kill ``a_4``; the surfaces would have done it anyway.)
* ``a_1 = \dfrac{p\,a^2b^2 + 2a_3\,a^2b^2\ln(b/a)}{b^2 - a^2}`` and ``a_2 = \ldots``: the Lamé solution of Section 7 **plus** a one-parameter family in ``a_3`` that leaves both curved surfaces traction-free. That family is the **pure bending of a curved beam** of Sadd §8.4.11, (8.4.60)–(8.4.61): its only resultant on a radial section is a moment.
* ``b_1 = -\dfrac{b_3}{2(a^2 + b^2)}``, ``b_2 = \dfrac{a^2b^2\,b_3}{2(a^2 + b^2)}``: the ``\sin\theta`` family collapses to a one-parameter family too, ``\sigma_r = -\dfrac{b_3}{a^2+b^2}\left(r + \dfrac{a^2b^2}{r^3} - \dfrac{a^2+b^2}{r}\right)\sin\theta``: the **curved beam loaded by a force at its end** of Sadd Exercise 8-34 (Timoshenko & Goodier §33). Its only resultant on the foot sections is a radial (here horizontal) force.

So the two free coefficients are the two **redundants** of a fixed arch: the foot moment and the horizontal reaction. The vertical reaction is not free — equilibrium fixes it. Let us confirm by integrating the tractions over the right foot. Its outward normal is ``-\mathbf{e}_\theta``, so the support exerts on the arch the traction ``(t_r, t_\theta) = (-\tau_{r\theta}, -\sigma_\theta)`` at ``\theta = 0``:

```math
H = -\int_R^{R+t}\tau_{r\theta}\,dr, \qquad V = -\int_R^{R+t}\sigma_\theta\,dr, \qquad M = -\int_R^{R+t}\sigma_\theta\,(r - r_m)\,dr
```

(``H`` positive outward, ``V`` positive upward, ``M`` positive counter-clockwise.)

Two technicalities before we compute them. The coefficients now contain ``\ln R`` and ``\ln(R+t)``, and `Symbolics`' fraction simplifier trips over sums of many fractions with such terms. `delog` replaces the two logarithms by plain symbols, `over_common_den` puts a sum of fractions over one denominator by hand, and `rsimp` simplifies the result and puts the logarithms back. `∫foot` is `∫r` across the foot with the same protection.
"""

# ╔═╡ d4000000-0000-0000-0000-000000000047
begin
	@variables LR LRt                                      # stand-ins for log R and log(R + t)
	delog(expr) = substitute(expr, Dict(log(R) => LR, log(R + t) => LRt))
	relog(expr) = substitute(expr, Dict(LR => log(R), LRt => log(R + t)))

	function fraction(e)                                   # e  ->  (N, D) with e = N/D, recursively, no gcd anywhere
		e = Symbolics.value(e)
		SymbolicUtils.iscall(e) || return (Num(e), Num(1))
		op, args = SymbolicUtils.operation(e), SymbolicUtils.arguments(e)
		if op === (+)
			parts = fraction.(args)
			udens = Num[d for d in unique(Symbolics.value.(last.(parts)))]
			D = prod(udens; init = Num(1))
			N = sum(parts[i][1]*prod(Num(d) for d in udens if !isequal(d, parts[i][2]); init = Num(1)) for i in eachindex(parts))
			return (N, D)
		elseif op === (*)
			parts = fraction.(args)
			return (prod(first.(parts); init = Num(1)), prod(last.(parts); init = Num(1)))
		elseif op === (/)
			(N1, D1), (N2, D2) = fraction(args[1]), fraction(args[2])
			return (N1*D2, D1*N2)
		elseif op === (^) && Symbolics.value(args[2]) isa Integer && Symbolics.value(args[2]) > 0
			(N1, D1) = fraction(args[1])
			return (N1^args[2], D1^args[2])
		else
			return (Num(e), Num(1))
		end
	end
	over_common_den(expr) = (ND = fraction(expand(expr)); (expand(ND[1]), ND[2]))
	function split_fractions(expr)                         # shallow version: only the top-level sum, enough for `simplify`
		e = Symbolics.value(expand(expr))
		terms = SymbolicUtils.iscall(e) && SymbolicUtils.operation(e) === (+) ? SymbolicUtils.arguments(e) : [e]
		parts = [SymbolicUtils.iscall(tm) && SymbolicUtils.operation(tm) === (/) ? Tuple(Num.(SymbolicUtils.arguments(tm))) : (Num(tm), Num(1)) for tm in terms]
		udens = Num[d for d in unique(Symbolics.value.(last.(parts)))]
		N = sum(expand(parts[i][1]*prod(Num(d) for d in udens if !isequal(d, parts[i][2]); init = Num(1))) for i in eachindex(parts))
		return expand(N), prod(udens; init = Num(1))
	end
	rsimp(expr)      = (ND = split_fractions(delog(expr)); relog(simplify(ND[1]/ND[2])))
	iszero_rat(expr) = (N = Symbolics.value(over_common_den(delog(expr))[1]); N isa Number && iszero(N))   # exact zero test, no simplification needed

	∫foot(expr)   = delog(∫r(delog(expr), R, R + t))     # ∫_R^{R+t} … dr, logarithms kept symbolic
	on_foot(expr) = at_angle(expr, 0, 0, 1)               # evaluate on the right foot θ = 0
	r_mid = R + t/2
end

# ╔═╡ d4000000-0000-0000-0000-000000000048
begin
	σr_1f  = substitute(σr_tf,  surface_map)
	σθ_1f  = substitute(σθ_tf,  surface_map)
	τrθ_1f = substitute(τrθ_tf, surface_map)
	H_foot_fx = rsimp(∫foot(-on_foot(τrθ_1f)))
	V_foot_fx = rsimp(∫foot(-on_foot(σθ_1f)))
	M_foot_fx = rsimp(∫foot(-on_foot(σθ_1f)*(r - r_mid)))
	(H = H_foot_fx, V = V_foot_fx, M = M_foot_fx)
end

# ╔═╡ d4000000-0000-0000-0000-000000000049
md"""
```math
V = p\,(R + t), \qquad
H = h_H\, b_3, \quad h_H = \frac{(a^2 + b^2)\ln(b/a) - (b^2 - a^2)}{a^2 + b^2}, \qquad
M = M_L + m_M\, a_3, \quad m_M = -\frac{(b^2 - a^2)^2 - 4a^2b^2\ln^2(b/a)}{2\,(b^2 - a^2)}
```

``V`` is the statics of Section 7 again, whatever ``a_3, b_3`` are. ``M_L = -\dfrac{p\,ab}{2}\left[1 - \dfrac{2ab\ln(b/a)}{b^2 - a^2}\right]`` is the small moment of the eccentric Lamé thrust found in Section 7 (now with the sign convention of a support reaction). ``m_M`` is Sadd's ``-N/[2(b^2 - a^2)]`` from (8.4.60), ``h_H`` is Timoshenko's ``N`` divided by ``a^2 + b^2``. The check cell confirms these three closed forms:
"""

# ╔═╡ d4000000-0000-0000-0000-00000000004a
begin
	a_, b_ = R, R + t
	hH  = ((a_^2 + b_^2)*log(b_/a_) - (b_^2 - a_^2))/(a_^2 + b_^2)
	mM  = -((b_^2 - a_^2)^2 - 4a_^2*b_^2*log(b_/a_)^2)/(2(b_^2 - a_^2))
	M_L = -p*a_*b_/2*(1 - 2a_*b_*log(b_/a_)/(b_^2 - a_^2))
	logsplit(expr) = substitute(expr, Dict(log(b_/a_) => log(R + t) - log(R)))     # ln(b/a) = ln b − ln a, so delog can see it
	resultant_check = iszero_rat.([H_foot_fx - logsplit(hH)*b3, M_foot_fx - logsplit(M_L) - logsplit(mM)*a3])
end

# ╔═╡ d4000000-0000-0000-0000-00000000004b
md"""
### 11.3 Displacements of the ``\sin\theta`` family

The ``a_3`` family is axisymmetric and its displacements are already in Section 5. For the new family ``\phi = f(r)\sin\theta`` the strains have the form ``\varepsilon_r = g_1(r)\sin\theta``, ``\varepsilon_\theta = g_2(r)\sin\theta``, ``\gamma_{r\theta} = g_3(r)\cos\theta``, and the integration of Section 5 goes through with one twist.

**Step 1.** ``u_r = U(r)\sin\theta + f(\theta)`` with ``U = \int g_1\,dr``.

**Step 2.** ``\partial u_\theta/\partial\theta = r\varepsilon_\theta - u_r`` gives ``u_\theta = (U - r g_2)\cos\theta - F(\theta) + g(r)``, ``F' = f``.

**Step 3.** In the shear strain the ``\cos\theta`` terms give ``[U' - rg_2'] \cos\theta + (f' + F)/r + g' - g/r = g_3\cos\theta``, i.e.

```math
\frac{f'(\theta) + F(\theta)}{r} + g'(r) - \frac{g(r)}{r} = h(r)\cos\theta, \qquad h = g_3 - g_1 + r g_2'
```

For the ``r^3`` and ``1/r`` terms compatibility makes ``h = 0`` and we recover the rigid-body motion of Section 5. For the ``r\ln r`` term ``h = C_0/r`` with ``C_0 = -4(1-\nu^2)\,b_3/E``, a constant that the ``\theta``-equation must absorb: ``F'' + F = C_0\cos\theta`` is *resonant* and its solution is ``F = \tfrac{1}{2}C_0\,\theta\sin\theta``, ``f = F' = \tfrac{1}{2}C_0(\sin\theta + \theta\cos\theta)``. The ``\theta\sin\theta`` and ``\theta\cos\theta`` terms are the same kind of multivaluedness as the ``r\theta`` term of Section 6: a full ring could not carry ``b_3 \neq 0`` (which is why Sadd's Exercise 8-34 is a curved *beam* problem), but our half ring can. Altogether

```math
u_r^{(b)} = U(r)\sin\theta + \tfrac{1}{2}C_0\left(\sin\theta + \theta\cos\theta\right), \qquad
u_\theta^{(b)} = \left[U(r) - r g_2(r)\right]\cos\theta - \tfrac{1}{2}C_0\,\theta\sin\theta
```

`Symbolics` does the bookkeeping: read ``g_1, g_2, g_3`` off the strains, integrate ``U`` with `∫r`, and check that ``h`` is indeed a constant over ``r``:
"""

# ╔═╡ d4000000-0000-0000-0000-00000000004c
begin
	ϕ_s = (b1*r^3 + b2/r + b3*r*log(r))*sin(θ)
	σr_s, σθ_s, τrθ_s = stresses_polar(ϕ_s)
	g1 = substitute(simplify((1 + ν)/E * ((1 - ν)*σr_s - ν*σθ_s)), Dict(sin(θ) => 1))
	g2 = substitute(simplify((1 + ν)/E * ((1 - ν)*σθ_s - ν*σr_s)), Dict(sin(θ) => 1))
	g3 = substitute(simplify(2(1 + ν)/E * τrθ_s), Dict(cos(θ) => 1))
	U_s = ∫r(g1)
	V_s = simplify(expand(U_s - r*g2))
	C0  = simplify(expand(r*(g3 - g1 + r*expand_derivatives(Dr(g2)))))
	(U = U_s, U_minus_rg2 = V_s, C0 = C0)
end

# ╔═╡ d4000000-0000-0000-0000-00000000004d
begin
	ur_general_fixed = ur_general + U_s*sin(θ) + C0/2*(sin(θ) + θ*cos(θ))
	uθ_general_fixed = uθ_general + V_s*cos(θ) - C0/2*θ*sin(θ)
	(ur_general_fixed, uθ_general_fixed)
end

# ╔═╡ d4000000-0000-0000-0000-00000000004e
md"Strain check of the full general field (``b_4 = 0``, as the surfaces require) — each residual should be zero:"

# ╔═╡ d4000000-0000-0000-0000-00000000004f
strain_check_general_fixed = let
	σr_, σθ_, τrθ_ = stresses_polar(ϕ_trial + ϕ_s)
	εr_  = (1 + ν)/E * ((1 - ν)*σr_ - ν*σθ_)
	εθ_  = (1 + ν)/E * ((1 - ν)*σθ_ - ν*σr_)
	γrθ_ = 2(1 + ν)/E * τrθ_
	simplify.(expand.([
		expand_derivatives(Dr(ur_general_fixed)) - εr_,
		expand_derivatives((ur_general_fixed + Dθ(uθ_general_fixed))/r) - εθ_,
		expand_derivatives(Dθ(ur_general_fixed)/r + Dr(uθ_general_fixed) - uθ_general_fixed/r) - γrθ_,
	]))
end

# ╔═╡ d4000000-0000-0000-0000-000000000050
md"""
### 11.4 What the fixed feet say — the force method in disguise

Substitute the surface results into the displacement field. Five unknowns remain: ``a_3, b_3`` and the rigid-body constants ``u_{x0}, u_{y0}, \omega_0``.

**Symmetry first.** ``u_\theta(r, \pi/2) = 0`` for all ``r`` is linear in ``r``, so its two coefficients give ``u_{x0}`` and ``\omega_0`` directly:
"""

# ╔═╡ d4000000-0000-0000-0000-000000000051
begin
	ur_1f = substitute(ur_general_fixed, surface_map)
	uθ_1f = substitute(uθ_general_fixed, surface_map)
	crown = expand(delog(at_angle(uθ_1f, Π/2, 1, 0)))                           # c₀ + c₁ r
	c₀ = substitute(crown, Dict(r => 0))
	c₁ = expand(substitute(crown, Dict(r => 1)) - c₀)
	sym_map = Dict(ux0 => expand(ux0 + c₀), ω0 => expand(ω0 - c₁))              # c₀ = −ux0 + …,  c₁ = ω0 + …
end

# ╔═╡ d4000000-0000-0000-0000-000000000052
md"""
```math
u_{x0} = \frac{\pi(1-\nu^2)}{E}\,b_3, \qquad \omega_0 = -\frac{2\pi(1-\nu^2)}{E}\,a_3
```

These are the multivalued terms evaluated at the crown: the ``\theta\cos\theta`` displacement of the ``b_3`` field and the ``r\theta`` displacement of the ``a_3`` field must be cancelled there by a rigid translation and a rigid rotation. Evaluated back at the right foot (``\theta = 0``), where the multivalued terms themselves vanish, they are what the feet feel.

**The three foot conditions.** Apply the functionals ``\Delta, \Lambda, \Theta`` to the right foot. ``\Lambda`` contains ``u_{y0}`` and nothing else new — it only pins the vertical rigid translation. ``\Delta`` and ``\Theta`` are linear in ``a_3`` and ``b_3``:

```math
\Delta = \Delta_L + \Delta_a\,a_3 + \Delta_b\,b_3, \qquad \Theta = \Theta_L + \Theta_a\,a_3 + \Theta_b\,b_3
```
"""

# ╔═╡ d4000000-0000-0000-0000-000000000053
begin
	ur_sym = substitute(ur_1f, sym_map)
	uθ_sym = substitute(uθ_1f, sym_map)
	slide(ur_)    = ∫foot(on_foot(ur_))/t
	lift(uθ_)     = ∫foot(on_foot(uθ_))/t
	rotation(uθ_) = 12*∫foot(on_foot(uθ_)*(r - r_mid))/t^3
	Δ, Λ, Θ = slide(ur_sym), lift(uθ_sym), rotation(uθ_sym)
	Δ_L = rsimp(substitute(Δ, Dict(a3 => 0, b3 => 0)))
	Δ_a = rsimp(Symbolics.derivative(Δ, a3))
	Δ_b = rsimp(Symbolics.derivative(Δ, b3))
	Θ_L = rsimp(substitute(Θ, Dict(a3 => 0, b3 => 0)))
	Θ_a = rsimp(Symbolics.derivative(Θ, a3))
	Θ_b = relog(Symbolics.derivative(Θ, b3))        # the polynomial gcd inside `simplify` fails on this one; keep it as it comes
	(Δ_L = Δ_L, Δ_b = Δ_b, Θ_L = Θ_L, Θ_a = Θ_a, Δ_a = Δ_a)
end

# ╔═╡ d4000000-0000-0000-0000-000000000054
md"""
* ``\Theta_L = 0`` and ``\Delta_L = -\dfrac{(1+\nu)K}{E\,t}\left[(1-2\nu)\dfrac{b^2 - a^2}{2} + a^2\ln\dfrac{b}{a}\right] < 0``: the Lamé field (the smooth-support solution) does not rotate the foot, and slides it *inward* by the mean of the ``u_r`` of Section 8. This is the displacement the fixed support must undo.
* ``\Delta_b = \pi(1-\nu^2)/E`` and ``\Theta_a = -2\pi(1-\nu^2)/E``: the end-force field slides the foot *uniformly* (its ``u_r(r,0)`` is just ``u_{x0}``) and the pure-bending field rotates it *rigidly* (its ``u_\theta(r,0)`` is just ``u_{y0} + \omega_0 r``). Both come straight from the multivalued terms.
* ``\Delta_a`` and ``\Theta_b`` are the cross terms: the slide caused by bending and the rotation caused by the end force.

Divide by the resultants ``m_M`` and ``h_H`` to turn coefficients per unit ``a_3, b_3`` into **flexibilities** per unit moment and per unit force,

```math
f_{\Delta H} = \frac{\Delta_b}{h_H} = \frac{\pi(1-\nu^2)(a^2+b^2)}{E\,N_H}, \qquad
f_{\Theta M} = \frac{\Theta_a}{m_M} = \frac{4\pi(1-\nu^2)(b^2-a^2)}{E\,N}, \qquad
f_{\Delta M} = \frac{\Delta_a}{m_M}, \qquad f_{\Theta H} = \frac{\Theta_b}{h_H}
```

with ``N = (b^2-a^2)^2 - 4a^2b^2\ln^2(b/a)`` (Sadd 8.4.60) and ``N_H = (a^2+b^2)\ln(b/a) - (b^2-a^2)``. **Maxwell's reciprocal theorem** says the slide per unit moment should equal the rotation per unit force, ``f_{\Delta M} = f_{\Theta H}``. The two were computed from completely different integrals (one of ``u_r`` of the axisymmetric field, one of ``u_\theta`` of the ``\sin\theta`` field), so comparing them is a genuine check. Reciprocity is *exact* only for work-conjugate pairs, and the plain mean of ``u_r`` across the foot is conjugate to a *uniform* horizontal traction, not to the actual (non-uniform) shear traction of the end-force field; the same goes for the mean rotation and the bending stress. So we expect agreement to ``O(t^2/r_m^2)``, not to all digits — and that is what the ratio below shows, for the parameters of Section 10:
"""

# ╔═╡ d4000000-0000-0000-0000-000000000055
begin
	f_ΔH = rsimp(Δ_b/logsplit(hH))
	f_ΘM = rsimp(Θ_a/logsplit(mM))
	f_ΔM = rsimp(Δ_a/logsplit(mM))
	f_ΘH = Θ_b/logsplit(hH)
	num_prm(expr) = build_function(expr, p, R, t, E, ν, Π; expression = Val(false))(prm.p, prm.R, prm.t, prm.E, prm.ν, π)
	maxwell_ratio = num_prm(f_ΔM)/num_prm(f_ΘH)
	(f_ΔH = f_ΔH, f_ΘM = f_ΘM, f_ΔM = f_ΔM, maxwell_ratio = maxwell_ratio)
end

# ╔═╡ d4000000-0000-0000-0000-000000000056
md"""
With ``M_b = m_M a_3`` the moment of the bending field and ``H = h_H b_3``, the two fixed-foot conditions read

```math
\begin{aligned}
\Delta_L + f_{\Delta M}\,M_b + f_{\Delta H}\,H &= 0 \\
f_{\Theta M}\,M_b + f_{\Theta H}\,H &= 0
\end{aligned}
\qquad\Longrightarrow\qquad
H = -\frac{\Delta_L\, f_{\Theta M}}{f_{\Delta H}\,f_{\Theta M} - f_{\Delta M}^2}, \qquad
M_b = -\frac{f_{\Delta M}}{f_{\Theta M}}\,H
```

These are the **compatibility equations of the force method** for a two-fold indeterminate arch: the smooth-support solution is the primary structure, ``\Delta_L`` is its incompatible displacement, and the two redundants ``H, M_b`` are whatever it takes to close the gap. The only difference from the structural-analysis version is that the flexibilities are *exact* — no beam hypothesis, no neglect of shear and axial deformation. The total foot moment is ``M = M_L + M_b``.
"""

# ╔═╡ d4000000-0000-0000-0000-000000000057
begin
	H_fixed  = -Δ_L*f_ΘM/(f_ΔH*f_ΘM - f_ΔM^2)
	Mb_fixed = -(f_ΔM/f_ΘM)*H_fixed
	ab_map   = Dict(a3 => Mb_fixed/logsplit(mM), b3 => H_fixed/logsplit(hH))
	uy0_fixed = relog(substitute(-substitute(Λ, Dict(uy0 => 0)), ab_map))            # Λ = uy0 + … = 0
	fixed_map = merge(ab_map, Dict(uy0 => uy0_fixed), Dict(k => substitute(v, ab_map) for (k, v) in sym_map))
	coeff_map_fixed = merge(Dict(k => substitute(v, fixed_map) for (k, v) in surface_map), fixed_map)
	(H = H_fixed, M_b = Mb_fixed)
end

# ╔═╡ d4000000-0000-0000-0000-000000000058
md"""
### 11.5 The thin-arch check: Castigliano's flexibilities

Strength of materials would set this problem up the same way, with the half-arch as a quarter-circle cantilever clamped at the crown (by symmetry the crown neither rotates nor moves sideways) and loaded at the foot by ``H`` and ``M_b``. Its bending moment is ``M(\theta) = M_b + H r_m\sin\theta`` and Castigliano's theorem with ``U = \int_0^{\pi/2} M(\theta)^2\, r_m\,d\theta/(2EI)`` gives

```math
f_{\Theta M}^{\text{thin}} = \frac{\pi r_m}{2EI}, \qquad
f_{\Theta H}^{\text{thin}} = f_{\Delta M}^{\text{thin}} = \frac{r_m^2}{EI}, \qquad
f_{\Delta H}^{\text{thin}} = \frac{\pi r_m^3}{4EI}, \qquad
EI = \frac{E\,t^3}{12(1-\nu^2)}
```

(plane strain stiffens the beam by ``1/(1-\nu^2)``). Then ``H^{\text{thin}} = -\Delta_L\,\dfrac{4\pi EI}{(\pi^2 - 8)\,r_m^3}`` and ``M_b^{\text{thin}} = -\dfrac{2r_m}{\pi}H^{\text{thin}}``: the textbook fixed semicircular arch with a prescribed spread of its feet. Expanding the exact flexibilities for ``t \ll R`` (``N \to \tfrac{4}{3}t^4``, ``N_H \to \tfrac{2}{3}t^3/r_m``) reproduces exactly these expressions, so the elasticity solution contains the curved-beam force method as its thin limit. The numerical example compares the two.

### 11.6 The solution and its checks

Substitute everything back. The stresses are the Lamé field plus ``M_b`` times the pure-bending field plus ``H`` times the end-force field; the displacements likewise. Verify, exactly: tractions on both curved surfaces for all ``\theta``, the three weak conditions on **both** feet (the left foot was never imposed — symmetry should take care of it), and the resultants:
"""

# ╔═╡ d4000000-0000-0000-0000-000000000059
begin
	σr_fx  = substitute(σr_tf,  coeff_map_fixed)
	σθ_fx  = substitute(σθ_tf,  coeff_map_fixed)
	τrθ_fx = substitute(τrθ_tf, coeff_map_fixed)
	σz_fx  = ν*(σr_fx + σθ_fx)
	ur_fx  = substitute(ur_general_fixed, coeff_map_fixed)
	uθ_fx  = substitute(uθ_general_fixed, coeff_map_fixed)
	left(expr) = at_angle(expr, Π, 0, -1)
	fixed_checks = (
		# curved surfaces, for all θ, whatever a3 and b3 are
		surfaces = iszero_rat.([trig(substitute(σr_1f, Dict(r => R))), trig(substitute(σr_1f, Dict(r => R + t))) + p,
		                        trig(substitute(τrθ_1f, Dict(r => R))), trig(substitute(τrθ_1f, Dict(r => R + t)))]),
		# mirror symmetry: the left foot slides like the right one and lifts and rotates the opposite way, whatever a3, b3, uy0 are,
		# so the three conditions imposed on the right foot hold on the left foot too
		left_foot_mirror = iszero_rat.([∫foot(left(ur_sym)) - t*Δ, ∫foot(left(uθ_sym)) + t*Λ, 12∫foot(left(uθ_sym)*(r - r_mid)) + t^3*Θ]),
	)
end

# ╔═╡ d4000000-0000-0000-0000-00000000005a
md"""
### 11.7 Numerical example, fixed feet

Same concrete lining as Section 10. `Π` is now ``\pi``.
"""

# ╔═╡ d4000000-0000-0000-0000-00000000005b
begin
	args_fx = (r, θ, p, R, t, E, ν, Π)
	prm_vals = (prm.p, prm.R, prm.t, prm.E, prm.ν, π)
	fx(expr) = (f = build_function(expr, args_fx...; expression = Val(false)); (ri, θi) -> f(ri, θi, prm_vals...))
	Sθ_fx, Sr_fx, Trθ_fx, Ur_fx, Uθ_fx = fx.((σθ_fx, σr_fx, τrθ_fx, ur_fx, uθ_fx))
	scalar(expr) = build_function(expr, p, R, t, E, ν, Π; expression = Val(false))(prm_vals...)
	H_num, Mb_num, ML_num = scalar(H_fixed), scalar(Mb_fixed), scalar(M_L)
	ΔL_num = scalar(Δ_L)
	EI = prm.E*prm.t^3/(12(1 - prm.ν^2)); rm_ = prm.R + prm.t/2
	thin = (fΘM = π*rm_/(2EI), fΘH = rm_^2/EI, fΔH = π*rm_^3/(4EI))
	H_thin  = -ΔL_num*4π*EI/((π^2 - 8)*rm_^3)
	Mb_thin = -2rm_/π*H_thin
	flex_max_dev = maximum(abs.([scalar(f_ΘM)/thin.fΘM, scalar(f_ΘH)/thin.fΘH, scalar(f_ΔH)/thin.fΔH] .- 1))
	# the final field, checked numerically across both feet with the trapezoid rule
	rs_foot = range(prm.R, prm.R + prm.t, length = 2001)
	∫trap(f) = (h = step(rs_foot); h*(sum(f, rs_foot) - (f(rs_foot[1]) + f(rs_foot[end]))/2))
	foot_checks = (
		right_foot = (slide = ∫trap(x -> Ur_fx(x, 0.0))/prm.t, lift = ∫trap(x -> Uθ_fx(x, 0.0))/prm.t,
		              rotation = 12∫trap(x -> Uθ_fx(x, 0.0)*(x - rm_))/prm.t^3),
		left_foot  = (slide = ∫trap(x -> Ur_fx(x, Float64(π)))/prm.t, lift = ∫trap(x -> Uθ_fx(x, Float64(π)))/prm.t,
		              rotation = 12∫trap(x -> Uθ_fx(x, Float64(π))*(x - rm_))/prm.t^3),
		resultants = (H = ∫trap(x -> -Trθ_fx(x, 0.0)) - H_num, V = ∫trap(x -> -Sθ_fx(x, 0.0)) - prm.p*(prm.R + prm.t),
		              M = ∫trap(x -> -Sθ_fx(x, 0.0)*(x - rm_)) - ML_num - Mb_num),
	)
	(H = H_num, M_b = Mb_num, M_L = ML_num, Δ_L = ΔL_num, H_thin = H_thin, Mb_thin = Mb_thin, checks = foot_checks)
end

# ╔═╡ d4000000-0000-0000-0000-00000000005c
fixed_comparison = let
	ri, ro = prm.R, prm.R + prm.t
	fmt(v) = string(round(v, sigdigits = 4) + 0.0)
	rows = [
		("Lamé foot slide Δ_L  [mm]",                     ΔL_num,                      "—"),
		("flexibility fΘM  [1/(N·mm)]",                  scalar(f_ΘM),                fmt(thin.fΘM)),
		("flexibility fΘH ≈ fΔM  [1/N]",                scalar(f_ΘH),                fmt(thin.fΘH)),
		("flexibility fΔH  [mm/N]",                      scalar(f_ΔH),                fmt(thin.fΔH)),
		("horizontal reaction H  [N/mm]",                 H_num,                       fmt(H_thin)),
		("bending-field foot moment M_b  [N·mm/mm]",      Mb_num,                      fmt(Mb_thin)),
		("Lamé eccentricity moment M_L  [N·mm/mm]",       ML_num,                      "—"),
		("vertical reaction V = p(R+t)  [N/mm]",          prm.p*ro,                    fmt(prm.p*ro)),
		("σθ foot, inner face  [MPa]",                    Sθ_fx(ri, 0.0),              fmt(Sθ(ri))),
		("σθ foot, outer face  [MPa]",                    Sθ_fx(ro, 0.0),              fmt(Sθ(ro))),
		("σθ crown, inner face  [MPa]",                   Sθ_fx(ri, π/2),              fmt(Sθ(ri))),
		("σθ crown, outer face  [MPa]",                   Sθ_fx(ro, π/2),              fmt(Sθ(ro))),
		("τrθ foot, mid-thickness  [MPa]",                Trθ_fx(rm_, 0.0),            fmt(0.0)),
		("uᵣ crown, inner face  [mm]",                    Ur_fx(ri, π/2),              fmt(Ur(ri))),
		("convergence of the inner diameter  [mm]",       2abs(Ur_fx(ri, π/2)),        fmt(2abs(Ur(ri)))),
		("largest uᵣ left on the foot (absolute value)  [mm]",               maximum(abs(Ur_fx(x, 0.0)) for x in range(ri, ro, length = 101)), fmt(abs(Ur(ri)))),
	]
	hdr = "| quantity | fixed feet (elasticity) | thin arch / smooth feet |\n|---|---|---|\n"
	Markdown.parse(hdr * join(["| $(q) | $(fmt(a)) | $(b) |" for (q, a, b) in rows], "\n"))
end

# ╔═╡ d4000000-0000-0000-0000-00000000005d
md"""
The third column holds the thin-arch flexibilities and reactions in the first block, and the smooth-support (Lamé) results in the second.

* The exact flexibilities agree with Castigliano's to within $(round(100flex_max_dev, digits = 2)) % for ``t/r_m = `` $(round(prm.t/rm_, digits = 3)), and ``H`` and ``M_b`` to within $(round(100*abs(H_thin/H_num - 1), digits = 1)) %. The curved-beam force method is an excellent approximation here; the elasticity solution tells us *how* good.
* The horizontal reaction ``H`` = $(round(H_num, digits = 2)) N/mm is only $(round(100*H_num/(prm.p*(prm.R + prm.t)), digits = 2)) % of the vertical one: holding the feet against a slide of $(round(abs(ΔL_num), digits = 3)) mm costs very little force. But the *bending* that force and the foot moment produce is not small: ``M_b`` = $(round(Mb_num, digits = 0)) N·mm/mm at the feet and ``M_b + H r_m`` = $(round(Mb_num + H_num*rm_, digits = 0)) N·mm/mm of the opposite sign at the crown, on top of the eccentricity moment ``M_L`` = $(round(ML_num, digits = 0)) N·mm/mm that the Lamé thrust already carried. The extreme-fibre bending stress ``6M_b/t^2`` ≈ $(round(6abs(Mb_num)/prm.t^2, digits = 2)) MPa is $(round(100*6abs(Mb_num)/prm.t^2/abs(Sθ(prm.R)), digits = 0)) % of the Lamé hoop stress on the inner face: the hoop stress at the feet rises from $(round(Sθ(prm.R), digits = 2)) to $(round(Sθ_fx(prm.R, 0.0), digits = 2)) MPa on the inner face, and at the crown it *drops* to $(round(Sθ_fx(prm.R, π/2), digits = 2)) MPa there while the outer face takes over. The support idealisation matters at the 10–20 % level — Saint-Venant tells us where the two solutions differ *beyond* that (within ``t`` of the feet), not that they agree.
* Fixing the feet does **not** stiffen the crown. The convergence nearly doubles, from $(round(2abs(Ur(prm.R)), digits = 2)) to $(round(2abs(Ur_fx(prm.R, π/2)), digits = 2)) mm: an arch whose feet are held apart while it wants to shrink must sag. (Castigliano on the thin arch gives the same answer: the foot of the quarter-ring cantilever rises by ``0.137\,H r_m^3/EI`` ≈ $(round(0.1366*H_num*rm_^3/EI, digits = 2)) mm relative to the crown.)
* The last row is what Saint-Venant leaves behind: the foot is fixed *on average*, but pointwise ``u_r(r, 0)`` still varies across it by a few hundredths of a millimetre, compared with the $(round(abs(ΔL_num), digits = 2)) mm slide of the smooth support. That residual, and the clamped-corner singularity of the exact solution, live within about ``t`` of the feet.
"""

# ╔═╡ d4000000-0000-0000-0000-00000000005e
fig_fixed = let
	θs = range(0, π, length = 181)
	ri, ro = prm.R, prm.R + prm.t
	fig = Figure(size = (1100, 420))

	ax1 = Axis(fig[1, 1], title = "Hoop stress around the arch", xlabel = "θ [°]", ylabel = "σθ [MPa]")
	lines!(ax1, rad2deg.(θs), Sθ_fx.(ri, θs), color = :firebrick, linewidth = 3, label = "inner face, fixed feet")
	lines!(ax1, rad2deg.(θs), Sθ_fx.(ro, θs), color = :steelblue, linewidth = 3, label = "outer face, fixed feet")
	hlines!(ax1, [Sθ(ri)], color = :firebrick, linestyle = :dash, linewidth = 2, label = "inner face, smooth feet")
	hlines!(ax1, [Sθ(ro)], color = :steelblue, linestyle = :dash, linewidth = 2, label = "outer face, smooth feet")
	axislegend(ax1, position = :rc)

	ax2 = Axis(fig[1, 2], title = "What is left on the right foot (smooth-support slide: $(round(ΔL_num, digits = 3)) mm)",
	           xlabel = "r [mm]", ylabel = "displacement [mm]")
	rs = range(ri, ro, length = 101)
	lines!(ax2, rs, Ur_fx.(rs, 0.0), color = :firebrick, linewidth = 3, label = "uᵣ(r, 0), fixed in the mean")
	lines!(ax2, rs, Uθ_fx.(rs, 0.0), color = :steelblue, linewidth = 3, label = "uθ(r, 0), fixed in the mean")
	hlines!(ax2, [0.0], color = (:gray, 0.5))
	axislegend(ax2, position = :rb)
	fig
end

# ╔═╡ d4000000-0000-0000-0000-00000000005f
md"""
Left: with fixed feet the hoop stress varies around the arch — the foot moment adds compression to the inner face at the feet and the crown moment, of the opposite sign, relieves it there and loads the outer face instead. Right: the price of the Saint-Venant relaxation. On the foot the mean of ``u_r`` and the mean and first moment of ``u_\theta`` are zero, but neither displacement is zero pointwise — compare the scale with the slide of the smooth support given in the title.

Deformed shapes, same exaggeration as Section 10. The fixed-foot arch is held at its feet and bends; the smooth-foot arch shrinks uniformly and slides.
"""

# ╔═╡ d4000000-0000-0000-0000-000000000060
fig_fixed_shape = let
	scale = 1500
	fig = Figure(size = (1100, 480))
	ax = Axis(fig[1, 1], title = "Deformed shape, displacements × $scale", xlabel = "x [mm]", ylabel = "y [mm]", aspect = DataAspect())
	θs = range(0, π, length = 181)
	ri, ro = prm.R, prm.R + prm.t
	# polar displacement -> cartesian: ux = ur cosθ − uθ sinθ,  uy = ur sinθ + uθ cosθ
	function arc_fx(rad, s)
		[Point2f(rad*cos(θi) + s*(Ur_fx(rad, θi)*cos(θi) - Uθ_fx(rad, θi)*sin(θi)),
		         rad*sin(θi) + s*(Ur_fx(rad, θi)*sin(θi) + Uθ_fx(rad, θi)*cos(θi))) for θi in θs]
	end
	foot_fx(s, θi) = [Point2f(x*cos(θi) + s*(Ur_fx(x, θi)*cos(θi) - Uθ_fx(x, θi)*sin(θi)),
	                          x*sin(θi) + s*(Ur_fx(x, θi)*sin(θi) + Uθ_fx(x, θi)*cos(θi))) for x in range(ri, ro, length = 11)]
	arc_sm(rad, s) = [Point2f((rad + s*Ur(rad))*cos(θi), (rad + s*Ur(rad))*sin(θi)) for θi in θs]
	for rad in (ri, ro)
		lines!(ax, arc_sm(rad, 0), color = (:gray, 0.8), linewidth = 2, label = rad == ri ? "undeformed" : nothing)
		lines!(ax, arc_sm(rad, scale), color = :steelblue, linewidth = 2, linestyle = :dash, label = rad == ri ? "smooth feet" : nothing)
		lines!(ax, arc_fx(rad, scale), color = :firebrick, linewidth = 3, label = rad == ri ? "fixed feet" : nothing)
	end
	for θi in (0.0, Float64(π))
		lines!(ax, foot_fx(0, θi), color = (:gray, 0.8), linewidth = 2)
		lines!(ax, foot_fx(scale, θi), color = :firebrick, linewidth = 3)
	end
	hlines!(ax, [0.0], color = :black, linewidth = 1)
	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ d4000000-0000-0000-0000-00000000003d
md"""
## 12. Summary

| quantity | symbolic (``K = \dfrac{p(R+t)^2}{t(2R+t)}``) | numeric |
|---|---|---|
| ``\phi`` | ``a_1 \ln r + a_2 r^2``, ``\;a_1 = K R^2``, ``\;a_2 = -K/2`` | ``K`` = $(round(K_num, digits = 3)) MPa |
| ``\sigma_r`` | ``-K\left(1 - R^2/r^2\right)`` | 0 at ``r_i``, ``-p`` = $(-prm.p) MPa at ``r_o`` |
| ``\sigma_\theta`` | ``-K\left(1 + R^2/r^2\right)`` | $(round(Sθ(prm.R), digits = 3)) MPa at ``r_i``, $(round(Sθ(prm.R + prm.t), digits = 3)) MPa at ``r_o`` |
| ``\tau_{r\theta}`` | ``0`` | 0 |
| ``\sigma_z`` | ``-2\nu K`` | $(round(Sz(prm.R), digits = 3)) MPa |
| hoop force per foot | ``N = \int_R^{R+t}\sigma_\theta\,dr = -p(R+t)`` | $(round(-prm.p*(prm.R + prm.t), digits = 1)) N/mm |
| ``u_r`` | ``-\dfrac{(1+\nu)K}{E}\left[(1-2\nu)\,r + \dfrac{R^2}{r}\right]`` | $(round(Ur(prm.R), digits = 4)) mm at ``r_i``, $(round(Ur(prm.R + prm.t), digits = 4)) mm at ``r_o`` |
| ``u_\theta``, ``u_z`` | ``0`` | 0 |
| thin ring ``\sigma_\theta`` | ``-p(R+t)/t`` (= exact average) | $(round(σθ_thin_num, digits = 3)) MPa |
| thin ring ``u_r`` | ``-(1-\nu^2)p(R+t)(R+t/2)/(Et)`` | $(round(ur_thin_num, digits = 4)) mm |

- The **polar Airy function** turns the circular geometry into a one-dimensional problem: of the Michell terms with ``\theta``-independent stresses, ``a_1\ln r + a_2 r^2 + a_3 r^2\ln r + a_4\theta``, the pressure fixes ``a_1, a_2`` and the **supports** fix ``a_3 = a_4 = 0`` — the ``r^2\ln r`` "dislocation" term would lift one foot relative to the other, the ``\theta`` term would shear the feet.
- Smooth supports at the feet are exactly the symmetry conditions of a **full ring** under uniform pressure, so the arch solution is Lamé's thick-walled cylinder (Sadd Example 8-6) in plane strain. What is relaxed is the friction at the feet and the non-uniformity of real ground pressure; Saint-Venant confines the first to the feet, the second needs the ``\cos 2\theta`` Michell terms.
- **Plane strain** leaves the stresses ``\sigma_r, \sigma_\theta`` unchanged from plane stress (traction problem, no body force), adds the uniform longitudinal stress ``\sigma_z = \nu(\sigma_r + \sigma_\theta)``, and stiffens the in-plane response: the displacement factor is ``(1+\nu)(1-2\nu)`` on the hydrostatic part instead of ``(1-\nu)``.
- The **thin-ring** hoop stress ``-p(R+t)/t`` is exactly the thickness-average of the true one; the true hoop stress is largest on the free inner face, by a factor ``2(R+t)^2/[(R+t)^2 + R^2]`` relative to the outer face, which ranges from 1 (thin lining) to 2 (hole in an infinite medium).
- **Fixed feet** (Section 11) make the problem a mixed one with no closed-form solution: a finite sum of Michell terms cannot clamp a foot pointwise, and the exact solution has clamped–free corner singularities. Fixing the feet **in the mean** (zero mean slide, lift and rotation — Saint-Venant, as at every beam end in Sadd Chapter 8) is solvable in closed form: the surfaces leave two free coefficients, ``a_3`` (Sadd's pure bending of a curved beam, §8.4.11) and ``b_3`` (the end-loaded curved beam, Exercise 8-34), which are the two redundants of a fixed arch, the foot moment ``M_b`` and the horizontal reaction ``H``. The foot conditions become the **force-method compatibility equations** with exact flexibilities whose thin limits are Castigliano's ``\pi r_m/2EI,\ r_m^2/EI,\ \pi r_m^3/4EI``.
- For the example lining ``H`` is only $(round(100*H_num/(prm.p*(prm.R + prm.t)), digits = 1)) % of ``V``, but the bending it induces changes the hoop stress by 10–20 % (``\sigma_\theta`` at the inner face: $(round(Sθ_fx(prm.R, 0.0), digits = 2)) MPa at the feet, $(round(Sθ_fx(prm.R, π/2), digits = 2)) MPa at the crown, against $(round(Sθ(prm.R), digits = 2)) MPa with smooth feet) and nearly doubles the convergence ($(round(2abs(Ur_fx(prm.R, π/2)), digits = 2)) vs $(round(2abs(Ur(prm.R)), digits = 2)) mm): an arch whose feet are held apart sags. The thin-arch force method reproduces ``H`` and ``M_b`` to within $(round(100*abs(H_thin/H_num - 1), digits = 1)) %.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
CairoMakie = "~0.15.15"
Latexify = "~0.16.12"
Symbolics = "~7.42.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.2"
manifest_format = "2.1"
project_hash = "be6ea7061fc7e9215d9dd51b691b4714ccbdf83e"

[[deps.ADTypes]]
deps = ["PrecompileTools"]
git-tree-sha1 = "629de23e1c16911b439dabd2303c08af9575b226"
registries = "General"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.24.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
registries = "General"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "e71ee7b4aa06b045259a7d6101e1cb45ad140bce"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.1"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
registries = "General"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "7063ad1083578215c7c4bf410368150abe8d5524"
registries = "General"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.45"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7c2c19b5a26e601634bf718490b89d59685f122e"
registries = "General"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.7.1"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
registries = "General"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
registries = "General"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
registries = "General"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
registries = "General"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "1aec1ff0dcaa83a484e7f72097562fa3102c6722"
registries = "General"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.30.2"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceAMDGPUExt = "AMDGPU"
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceFillArraysExt = "FillArrays"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceGPUArraysCoreTrackerExt = ["GPUArraysCore", "Tracker"]
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FillArrays = "1a297f60-69ca-5386-bcde-b61e274b549b"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "94eab0b3ccdcac361188cc661daf69d4433c1818"
registries = "General"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.2.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
registries = "General"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
registries = "General"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "8c290a1b223deaeea9aea44b235d24546da8eb98"
registries = "General"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.4.0"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
registries = "General"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
registries = "General"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
registries = "General"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
registries = "General"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
registries = "General"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
registries = "General"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "1cda0b7d5abfc95357dae18aca934d401f7869ad"
registries = "General"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.15.15"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "7b841680738c19948120f6e4cf8d0200518bd564"
registries = "General"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.8+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "12177ad6b3cad7fd50c8b3825ce24a99ad61c18f"
registries = "General"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZstd]]
deps = ["TranscodingStreams", "Zstd_jll"]
git-tree-sha1 = "da54a6cd93c54950c15adf1d336cfd7d71f51a56"
registries = "General"
uuid = "6b39b394-51ab-5f42-8807-6242bab2b4c2"
version = "0.8.7"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
registries = "General"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
registries = "General"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "61761f58648aa7217445f24f841839b78c712232"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.3"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
registries = "General"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Reexport"]
git-tree-sha1 = "291665b547f137df070e4dd83e432b5fee8cc4a0"
registries = "General"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.2"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
registries = "General"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
deps = ["PrecompileTools"]
git-tree-sha1 = "6c389fa857f6ca5a95474b52a52023fd77f24cb7"
registries = "General"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.14"

[[deps.CommonWorldInvalidations]]
deps = ["PrecompileTools"]
git-tree-sha1 = "7c4ea0adfb7238bf4678aa36325863fab6163f6f"
registries = "General"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.2.2"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
registries = "General"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
registries = "General"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
registries = "General"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "7bc84b769c1d384315e7b5c4ac03a6c303e6cf35"
registries = "General"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.8"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
registries = "General"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
registries = "General"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CoreMath]]
deps = ["CoreMath_jll"]
git-tree-sha1 = "8c0480f92b1b1796239156a1b9b1bfb1b39499b4"
registries = "General"
uuid = "b7a15901-be09-4a0e-87d2-2e66b0e09b5a"
version = "0.1.0"

[[deps.CoreMath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a692a4c1dc59a4b8bc0b6403876eb3250fde2bc3"
registries = "General"
uuid = "a38c48d9-6df1-5ac9-9223-b6ada3b5572b"
version = "0.1.0+0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
registries = "General"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
registries = "General"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
registries = "General"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "4ac548adcad90c1d5d677af13568a748af4c952b"
registries = "General"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.7"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "79a2aca180a85c690c58a020d47b426954b590f8"
registries = "General"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.16.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "Roots", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "a958ab3a40c755563f5e1405c0846cb0446bf19d"
registries = "General"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.131"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
registries = "General"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.DomainSets]]
deps = ["CompositeTypes", "FunctionMaps", "IntervalSets", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "6ac879a64381881dc24f96a15e13e80204439906"
registries = "General"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.8.3"
weakdeps = ["Makie", "Random", "RecipesBase"]

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"
    DomainSetsRandomExt = "Random"
    DomainSetsRecipesBaseExt = "RecipesBase"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.DynamicPolynomials]]
deps = ["LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "StarAlgebras"]
git-tree-sha1 = "b95c4325301d86a6a1b59790b71f7269a8702dd5"
registries = "General"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
registries = "General"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
registries = "General"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
registries = "General"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2bfb1e047e2ad0a5ca94365340bde8005d637568"
registries = "General"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.8.4+0"

[[deps.ExprTools]]
git-tree-sha1 = "d2e49e7efd29719d6f28b891b0e0e159daa9d2b4"
registries = "General"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.11"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
registries = "General"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d9d3cd382f2c999f684e6bfa4039fcc16c153b6f"
registries = "General"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "9.0.2+0"

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "65e55303b72f4a567a51b174dd2c47496efeb95a"
registries = "General"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6621fef488e496356c9c9625d0562c12a6070819"
registries = "General"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.20.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport"]
git-tree-sha1 = "a1b2fbfe98503f15b665ed45b3d149e5d8895e4c"
registries = "General"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.9.0"

    [deps.FilePaths.extensions]
    FilePathsGlobExt = "Glob"
    FilePathsURIParserExt = "URIParser"
    FilePathsURIsExt = "URIs"

    [deps.FilePaths.weakdeps]
    Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
    URIParser = "30578b45-9adc-5946-b283-645ec420af67"
    URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
registries = "General"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "086b5fbd032baf544678cc15b52b015e4c4aceb8"
registries = "General"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.17.1"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
registries = "General"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
registries = "General"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
registries = "General"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
registries = "General"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
registries = "General"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
registries = "General"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.FunctionMaps]]
deps = ["CompositeTypes", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "31bd99a57edf98990d1c21486032963955450e8d"
registries = "General"
uuid = "a85aefff-f8ca-4649-a888-c8e5398bc76c"
version = "0.1.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.Gamma]]
deps = ["LogExpFunctions"]
git-tree-sha1 = "becc397f7cfb06e343496ae6ffb04818a851da51"
registries = "General"
uuid = "a0844989-3bd2-4988-8bea-c9407ab0941b"
version = "1.2.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "ec46c5825710fa1a15d468acb2d93cc939a7a5fe"
registries = "General"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.13"

    [deps.GeometryBasics.extensions]
    ExtentsExt = "Extents"
    GeometryBasicsGeoInterfaceExt = "GeoInterface"
    IntervalSetsExt = "IntervalSets"

    [deps.GeometryBasics.weakdeps]
    Extents = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
registries = "General"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
registries = "General"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3efbbc027441271444dcd0c0a46f2d119dc4329"
registries = "General"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "6.1.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "090526e65de8f69648ac156daae153de8b56df62"
registries = "General"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.88.3+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
registries = "General"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "69ffb934a5c5b7e086a0b4fee3427db2556fba6e"
registries = "General"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.16+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Inflate", "LinearAlgebra", "Random", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "dbb2b976f605b220dfe139d6db9f3859680da09e"
registries = "General"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.15.0"
weakdeps = ["Distributed", "SharedArrays"]

    [deps.Graphs.extensions]
    GraphsSharedArraysExt = "SharedArrays"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "ef70da5e123a06a29e2d6ddff0f09985bc226491"
registries = "General"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "9d9531a9cb63a9edc33836414e82a07e81710de2"
registries = "General"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "100.14004.0+0"

[[deps.HypergeometricFunctions]]
deps = ["Gamma", "LinearAlgebra"]
git-tree-sha1 = "31bb6c92405c084617facc1d7ed9eb6c402d061e"
registries = "General"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.30"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
registries = "General"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
registries = "General"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
registries = "General"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "f0f005f997dfb8c5fe23920d99458a9619873893"
registries = "General"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.10"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
registries = "General"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dcc8d0cd653e55213df9b75ebc6fe4a8d3254c65"
registries = "General"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.2.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
registries = "General"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
registries = "General"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "c72458f1962faeb003bf23cbdb75164fe6280906"
registries = "General"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "48922d06068130f87e43edef52382e6a94305ae6"
registries = "General"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.3"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "CoreMath", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "1c531bf0f8a5c60a340926e058fd3f209b5eef5d"
registries = "General"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.12"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticIrrationalConstantsExt = "IrrationalConstants"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticMakieExt = "Makie"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    IrrationalConstants = "92d709cd-6900-40b7-9082-c6be49f344b6"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "0db6aea5b64caa1c47e9fdd27394f0296f05e8bf"
registries = "General"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.15"
weakdeps = ["Makie", "Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsMakieExt = "Makie"
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
registries = "General"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
registries = "General"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
registries = "General"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
registries = "General"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
registries = "General"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
registries = "General"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "4657a834b01ce00f2b94546095c4350617b9af2c"
registries = "General"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.10.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
registries = "General"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
registries = "General"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "037babc10853eeb8e585418922246cb97b8e5b74"
registries = "General"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.2.0+1"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTA", "Interpolations", "StatsBase"]
git-tree-sha1 = "9eda8292dd3268b3b7ec9df21bbfac24e177ec52"
registries = "General"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.12"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
registries = "General"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "39bca05343661c347aae0bca57a5994a0bf4f08d"
registries = "General"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.2.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e5b100780d4d30d63b4618d7930d48af409c1772"
registries = "General"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "23.1.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f88f3ccef05a6a72a0cf0ed417c8fd68530f4ab2"
registries = "General"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.1"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "df7566479bd64f20bd16b09960145e70160ffb3b"
registries = "General"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.12"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
registries = "General"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "1.0.0"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.22.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.7+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "OpenSSL_jll", "Zlib_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.104+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
registries = "General"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
registries = "General"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
registries = "General"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc3ad4faf30015a3e8094c9b5b7f19e85bdf2386"
registries = "General"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "aebd334d06cee9f24cea70bd19a39749daf73881"
registries = "General"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.3+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d620582b1f0cbe2c72dd1d5bd195a9ce73370ab1"
registries = "General"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.13.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "b85e2797b2409570e84c4de46238c0ed5f6476ae"
registries = "General"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "1.0.2"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
registries = "General"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "5f6f5d1b1fb7ff98c9a083bbfd4c661a9808e758"
registries = "General"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.15"

    [deps.Makie.extensions]
    MakieDynamicQuantitiesExt = "DynamicQuantities"

    [deps.Makie.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"

[[deps.MappedArrays]]
git-tree-sha1 = "0ee4497a4e80dbd29c058fcee6493f5219556f40"
registries = "General"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.3"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "aa1078778be5a8e5259ff04fbc3d258b3e78d464"
registries = "General"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.9"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
registries = "General"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
registries = "General"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "3f234cf937eb1d58f613440ecb022da784ad34d4"
registries = "General"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.13"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.8.13"

[[deps.MuladdMacro]]
deps = ["PrecompileTools"]
git-tree-sha1 = "283bf85d4a767481dd924dff0eee1735e95f449e"
registries = "General"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.7"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics", "StarAlgebras"]
git-tree-sha1 = "65f1160de342181533efa7061f86539d420ab81d"
registries = "General"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.20"
weakdeps = ["ChainRulesCore", "Latexify"]

    [deps.MultivariatePolynomials.extensions]
    LatexifyExt = "Latexify"
    MultivariatePolynomialsChainRulesCoreExt = "ChainRulesCore"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "7fc67976fca2c577ca96257061e4785463bb9fb8"
registries = "General"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.8.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "dbd2e8cd2c1c27f0b584f6661b4309609c5a685e"
registries = "General"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.4"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
registries = "General"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
registries = "General"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
registries = "General"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
registries = "General"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "38a93f17e431141c6470bb67a88952a7c4f0e928"
registries = "General"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.34+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
registries = "General"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e9ae72527609bb66882d38ad009ffd47ca5713fb"
registries = "General"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.4.16+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.8+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
registries = "General"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
registries = "General"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05f45c2e0de6259db764adbfd2f1dc6d3f8de13c"
registries = "General"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "2.0.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.46.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "123266c25174ef6c8d4718920abc206452cf8de6"
registries = "General"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.41"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "32b657a0d57c310a1a172bfc8c8cf68c5e674323"
registries = "General"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.5"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
registries = "General"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
registries = "General"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1912a9f1b9ca55005b03ba075f8e19993583e237"
registries = "General"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.58.2+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools"]
git-tree-sha1 = "663e8b48b789916221e0765393b289ca6c88f24e"
registries = "General"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "3.0.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "e4a6721aa89e62e5d4217c0b21bd714263779dda"
registries = "General"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.46.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
registries = "General"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Reexport", "Statistics"]
git-tree-sha1 = "f20e945b895d2009c6c28d8bbf40a5cd846f7c2f"
registries = "General"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.5.0"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
registries = "General"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
registries = "General"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "5005266de4bfe50e53ff44a5cb5c540b6e47a254"
registries = "General"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.6.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
registries = "General"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
registries = "General"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
registries = "General"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "472daaa816895cb7aee81658d4e7aec901fa1106"
registries = "General"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "5e8e8b0ab68215d7a2b14b9921a946fee794749e"
registries = "General"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.3"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["Base64", "Dates", "FileWatching", "InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
registries = "General"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
registries = "General"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.ReadOnlyArrays]]
git-tree-sha1 = "e6f7ddf48cf141cb312b078ca21cb2d29d0dc11d"
registries = "General"
uuid = "988b38a3-91fc-5605-94a2-ee2116b3bd83"
version = "0.2.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "0c0dcffc464cda2ff33fc8ffc1dd1bdd53cd98be"
registries = "General"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.4.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
registries = "General"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
registries = "General"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
registries = "General"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d40b2fe70437b01397d2a4d5b020008da4e7019"
registries = "General"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.2+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "971f04b3780c0da4edf230573fc1aadc2c6e649e"
registries = "General"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "3.0.10"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
registries = "General"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "PrecompileTools", "SHA", "Serialization"]
git-tree-sha1 = "9bb075002b40849582905394bb3031fb77aa102d"
registries = "General"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.27"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
registries = "General"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.SciMLPublic]]
deps = ["PrecompileTools"]
git-tree-sha1 = "74685afb51732a464fbce79a72708f7c4203ceb7"
registries = "General"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.3.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
registries = "General"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
registries = "General"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "57aa595158717ef165e6f5ab639fe2e3178c0a2b"
registries = "General"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SignedDistanceFields]]
deps = ["Statistics"]
git-tree-sha1 = "3949ad92e1c9d2ff0cd4a1317d5ecbba682f4b92"
registries = "General"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "7ddb0b49c109481b046972c0e4ab02b2127d6a75"
registries = "General"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2c79185e5261b159474903598c571ddd9f12ff7b"
registries = "General"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.6"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
registries = "General"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.13.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "429071b23f4c9a13fb6582f807cc2ef454082408"
registries = "General"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.9.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
registries = "General"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StarAlgebras]]
deps = ["LinearAlgebra", "MutableArithmetics", "SparseArrays"]
git-tree-sha1 = "235b1f9d287bbf34083b3d0829343a7942c0ad1c"
registries = "General"
uuid = "0c0c59c1-dc5f-42e9-9a8b-b5dc384a6cd1"
version = "0.3.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "39e70e0ab5d7f89833a62ab7c79df15d4fc417c1"
registries = "General"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.22"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
registries = "General"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e2b53ce13a53367e96601081e33d34746b571bad"
registries = "General"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.5"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
registries = "General"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "adb9da019510162e67a4493fc235c23203d8b09e"
registries = "General"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.13"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91a5737baed20ee31f3faea0e51f57461f6a689e"
registries = "General"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "2.2.1"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "ad8002667372439f2e3611cfd14097e03fa4bccd"
registries = "General"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.3"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "b814d5005d6a529d740ffe06f8a86396f6501138"
registries = "General"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.9.2"

    [deps.StructUtils.extensions]
    StructUtilsLazilyInitializedFieldsExt = ["LazilyInitializedFields"]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    LazilyInitializedFields = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.10.1+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "PrecompileTools", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "7eb6da9656581ac9fbffaef3ef7950a090a5002a"
registries = "General"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.55"

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

    [deps.SymbolicIndexingInterface.weakdeps]
    PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[[deps.SymbolicLimits]]
deps = ["PrecompileTools", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "3780b4d8aaa1db8996a67146ac7d4ac20606f56e"
registries = "General"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "1.2.1"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "EnumX", "ExproniconLite", "Graphs", "LinearAlgebra", "MacroTools", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "OrderedCollections", "PrecompileTools", "ReadOnlyArrays", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "WeakCacheSets"]
git-tree-sha1 = "5b6c4665d3c22dda9c991a9426e72e9fe02305e6"
registries = "General"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "4.49.0"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsChainRulesCoreExt = "ChainRulesCore"
    SymbolicUtilsDistributionsExt = "Distributions"
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "AbstractPlutoDingetjes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "IntervalSets", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "Preferences", "Primes", "RecipesBase", "RuntimeGeneratedFunctions", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "c039acf6cded9bf0cbcada1c1ef60ca7b309b22c"
registries = "General"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "7.42.0"

    [deps.Symbolics.extensions]
    SymbolicsD3TreesExt = "D3Trees"
    SymbolicsDistributionsExt = "Distributions"
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsHypergeometricFunctionsExt = "HypergeometricFunctions"
    SymbolicsLatexifyExt = ["Latexify", "LaTeXStrings"]
    SymbolicsNemoExt = "Nemo"
    SymbolicsStaticArraysExt = "StaticArrays"
    SymbolicsSymPyExt = "SymPy"
    SymbolicsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Symbolics.weakdeps]
    D3Trees = "e3df1716-f71e-5df9-9e2d-98e193103c45"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    HypergeometricFunctions = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
registries = "General"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "a94d9bdda1b7bed0046cea645639ab3f62196fac"
registries = "General"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.14.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
registries = "General"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
registries = "General"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
registries = "General"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TiffImages]]
deps = ["CodecZstd", "ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "9ca5f1f2d42f80df4b8c9f6ab5a64f438bbd9976"
registries = "General"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.9"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
registries = "General"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
registries = "General"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
registries = "General"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "1f0f9f401753701a7e4113b5056ca38d33875b55"
registries = "General"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.29.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.WeakCacheSets]]
git-tree-sha1 = "386050ae4353310d8ff9c228f83b1affca2f7f38"
registries = "General"
uuid = "d30d5f5c-d141-4870-aa07-aabb0f5fe7d5"
version = "0.1.0"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
registries = "General"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
registries = "General"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e52eca002a11c30a858185efdfb15311e1c7a6bf"
registries = "General"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
registries = "General"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
registries = "General"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
registries = "General"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
registries = "General"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
registries = "General"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
registries = "General"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "58972370b81423fc546c56a60ed1a009450177c3"
registries = "General"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.19.0+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
registries = "General"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
registries = "General"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.2+0"

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
registries = "General"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1210ba774d3427387d307bf1f416d699b7c39417"
registries = "General"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.15.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "cb007192783c56d8249db4cf0e3495001edfe414"
registries = "General"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.5+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "28e57478e8a160d346a19c28b3fffb9273bcc9c2"
registries = "General"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.134+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
registries = "General"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "32781be40fe86735af02eae0dd22754e4b5f779d"
registries = "General"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.59+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "e067c8bae65bb40866552296a2d5b4dc65b8928f"
registries = "General"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.80.702+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
registries = "General"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
registries = "General"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "52d3b9475133c3bc8c0a7f90f18bfc8cdb5a443a"
registries = "General"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.70.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.8.3+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
registries = "General"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
registries = "General"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

# ╔═╡ Cell order:
# ╠═d4000000-0000-0000-0000-000000000001
# ╟─d4000000-0000-0000-0000-000000000002
# ╟─d4000000-0000-0000-0000-000000000003
# ╟─d4000000-0000-0000-0000-000000000004
# ╟─d4000000-0000-0000-0000-000000000005
# ╠═d4000000-0000-0000-0000-000000000006
# ╠═d4000000-0000-0000-0000-000000000007
# ╠═d4000000-0000-0000-0000-000000000008
# ╠═d4000000-0000-0000-0000-000000000009
# ╟─d4000000-0000-0000-0000-00000000000a
# ╠═d4000000-0000-0000-0000-00000000000b
# ╠═d4000000-0000-0000-0000-00000000000c
# ╟─d4000000-0000-0000-0000-00000000000d
# ╠═d4000000-0000-0000-0000-00000000000e
# ╟─d4000000-0000-0000-0000-00000000000f
# ╠═d4000000-0000-0000-0000-000000000010
# ╟─d4000000-0000-0000-0000-000000000011
# ╠═d4000000-0000-0000-0000-000000000012
# ╟─d4000000-0000-0000-0000-000000000013
# ╠═d4000000-0000-0000-0000-000000000014
# ╟─d4000000-0000-0000-0000-000000000015
# ╠═d4000000-0000-0000-0000-000000000016
# ╟─d4000000-0000-0000-0000-000000000017
# ╠═d4000000-0000-0000-0000-000000000018
# ╟─d4000000-0000-0000-0000-000000000019
# ╠═d4000000-0000-0000-0000-00000000001a
# ╠═d4000000-0000-0000-0000-00000000001b
# ╟─d4000000-0000-0000-0000-00000000001c
# ╠═d4000000-0000-0000-0000-00000000001d
# ╠═d4000000-0000-0000-0000-00000000001e
# ╟─d4000000-0000-0000-0000-00000000001f
# ╟─d4000000-0000-0000-0000-000000000020
# ╠═d4000000-0000-0000-0000-000000000021
# ╠═d4000000-0000-0000-0000-000000000022
# ╟─d4000000-0000-0000-0000-000000000023
# ╠═d4000000-0000-0000-0000-000000000024
# ╠═d4000000-0000-0000-0000-000000000025
# ╟─d4000000-0000-0000-0000-000000000026
# ╠═d4000000-0000-0000-0000-000000000027
# ╟─d4000000-0000-0000-0000-000000000028
# ╠═d4000000-0000-0000-0000-000000000029
# ╟─d4000000-0000-0000-0000-00000000002a
# ╠═d4000000-0000-0000-0000-00000000002b
# ╠═d4000000-0000-0000-0000-00000000002c
# ╟─d4000000-0000-0000-0000-00000000002d
# ╠═d4000000-0000-0000-0000-00000000002e
# ╟─d4000000-0000-0000-0000-00000000002f
# ╠═d4000000-0000-0000-0000-000000000030
# ╟─d4000000-0000-0000-0000-000000000031
# ╠═d4000000-0000-0000-0000-000000000032
# ╟─d4000000-0000-0000-0000-00000000003e
# ╟─d4000000-0000-0000-0000-000000000033
# ╠═d4000000-0000-0000-0000-000000000034
# ╠═d4000000-0000-0000-0000-000000000035
# ╟─d4000000-0000-0000-0000-000000000036
# ╠═d4000000-0000-0000-0000-000000000037
# ╟─d4000000-0000-0000-0000-000000000038
# ╟─d4000000-0000-0000-0000-000000000039
# ╠═d4000000-0000-0000-0000-00000000003a
# ╟─d4000000-0000-0000-0000-00000000003b
# ╠═d4000000-0000-0000-0000-00000000003c
# ╟─d4000000-0000-0000-0000-00000000003f
# ╠═d4000000-0000-0000-0000-000000000040
# ╠═d4000000-0000-0000-0000-000000000041
# ╠═d4000000-0000-0000-0000-000000000042
# ╟─d4000000-0000-0000-0000-000000000043
# ╠═d4000000-0000-0000-0000-000000000044
# ╠═d4000000-0000-0000-0000-000000000045
# ╟─d4000000-0000-0000-0000-000000000046
# ╠═d4000000-0000-0000-0000-000000000047
# ╠═d4000000-0000-0000-0000-000000000048
# ╟─d4000000-0000-0000-0000-000000000049
# ╠═d4000000-0000-0000-0000-00000000004a
# ╟─d4000000-0000-0000-0000-00000000004b
# ╠═d4000000-0000-0000-0000-00000000004c
# ╠═d4000000-0000-0000-0000-00000000004d
# ╟─d4000000-0000-0000-0000-00000000004e
# ╠═d4000000-0000-0000-0000-00000000004f
# ╟─d4000000-0000-0000-0000-000000000050
# ╠═d4000000-0000-0000-0000-000000000051
# ╟─d4000000-0000-0000-0000-000000000052
# ╠═d4000000-0000-0000-0000-000000000053
# ╟─d4000000-0000-0000-0000-000000000054
# ╠═d4000000-0000-0000-0000-000000000055
# ╟─d4000000-0000-0000-0000-000000000056
# ╠═d4000000-0000-0000-0000-000000000057
# ╟─d4000000-0000-0000-0000-000000000058
# ╠═d4000000-0000-0000-0000-000000000059
# ╟─d4000000-0000-0000-0000-00000000005a
# ╠═d4000000-0000-0000-0000-00000000005b
# ╠═d4000000-0000-0000-0000-00000000005c
# ╟─d4000000-0000-0000-0000-00000000005d
# ╠═d4000000-0000-0000-0000-00000000005e
# ╟─d4000000-0000-0000-0000-00000000005f
# ╠═d4000000-0000-0000-0000-000000000060
# ╟─d4000000-0000-0000-0000-00000000003d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
