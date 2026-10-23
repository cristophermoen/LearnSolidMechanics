### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ e5000000-0000-0000-0000-000000000001
begin
	using Symbolics
	using CairoMakie
	using Latexify      # lets Pluto render Symbolics expressions as LaTeX
end

# ╔═╡ e5000000-0000-0000-0000-000000000002
md"""
# Thin Plate with a Circular Hole Under Uniaxial Tension — Airy Stress Function Solution

*Reference solution for the Neural Operator + FEA in-class project.*

A thin rectangular plate, ``L = 300`` mm long and ``W = 50`` mm wide, has a circular hole of diameter ``d = 25`` mm (radius ``a = 12.5`` mm) at its centre. Its two ends are pulled by a uniform tension ``T`` (force per unit area, i.e. end force divided by ``W`` times the thickness), the long edges and the hole are free, and the material is linear elastic with ``E, \nu``. The plate is thin, so this is **plane stress** (Sadd §7.2).

The hole is circular, so — as for the tunnel of Lecture 13 — the natural tool is the **polar-coordinate Airy stress function** and the **Michell solution** (Sadd §7.6, §8.3). The closed-form answer is the celebrated **Kirsch (1898) solution** for a stress-free hole in an *infinite* plate, Sadd Example 8-7, with its stress concentration factor of 3. The plate in the figure is not infinite: its width is only four hole radii. Section 9 measures exactly what the infinite-plate solution gets wrong on this strip, and that gap is one of the things the finite element model and the neural operator of the project are there to capture.

References: Sadd, *Elasticity*, §7.2 (plane stress), §7.6 (polar formulation), §8.3 (Michell solution), §8.4.1–8.4.2 (axisymmetric fields), **Example 8-7 and (8.4.11)–(8.4.17)** (stress-free hole under far-field tension), Figure 8-13. Finite-width correction: Howland, *Phil. Trans. R. Soc. A* 229 (1930); Pilkey & Pilkey, *Peterson's Stress Concentration Factors*, Chart 4.1.

**Coordinates.** Origin at the centre of the hole; ``x`` along the load (the long direction), ``y`` across the width; polar angle ``\theta`` measured from the ``+x`` axis. The plate occupies ``|x| \le L/2``, ``|y| \le W/2``, ``r \ge a``. The radial displacement ``u_r`` is positive outward and ``u_\theta`` positive counter-clockwise.

Steps:
1. Decide what in the figure can be solved in closed form and what cannot (the idealisation).
2. Write the governing equations in polar coordinates for plane stress.
3. State the boundary conditions, including the far-field tension in polar form.
4. Pick the terms of the Michell solution that the problem needs and impose boundedness at infinity.
5. Enforce the far-field and hole conditions symbolically with `Symbolics.jl`; recover the stress field and the stress concentration factor.
6. Integrate the strains to get the displacements; find the deformed shape of the hole.
7. Convert everything to Cartesian components (what a finite element code reports).
8. Work the numerical example with the dimensions in the figure; plot stresses, deformations and the residual tractions the infinite-plate solution leaves on the strip edges; quote the finite-width stress concentration factor.
"""

# ╔═╡ e5000000-0000-0000-0000-000000000003
md"""
## 1. Idealisation — from the figure to a solvable problem

The figure has three lengths: the hole radius ``a = 12.5`` mm, the half-width ``W/2 = 25`` mm ``= 2a`` and the half-length ``L/2 = 150`` mm ``= 12a``. An elasticity solution with a hole in a *rectangle* has no closed form: the hole wants polar coordinates, the edges want Cartesian ones, and no finite sum of Michell terms satisfies both sets of boundary conditions at once. What does have a closed form is the hole in an **infinite plate** under a uniform far-field tension ``T``. Whether that is a good model of the strip is a question of how fast the disturbance caused by the hole dies out.

We will find that the disturbance decays like ``(a/r)^2`` (and ``(a/r)^4``). At the **ends** ``x = \pm L/2 = \pm 12a`` that is ``1/144 < 1\,\%``: the ends are at infinity for all practical purposes, and the uniform traction ``T`` applied there *is* the far-field tension. At the **long edges** ``y = \pm W/2 = \pm 2a`` it is ``1/4``: the edges are not at infinity. The infinite-plate solution leaves non-zero normal and shear tractions on ``y = \pm W/2`` that the real free edges cannot supply, and Section 9 plots them. Removing them needs either Howland's (1930) infinite series or a numerical method — in the project, the finite element model.

So the plan is: solve the infinite plate **exactly** (Sections 2–8), then quantify its error on this strip and quote the accepted finite-width correction (Section 9). The infinite-plate solution is still the right starting point: it carries the physics (the ``3T`` hoop stress, the compressive ``-T`` at the sides of the hole, the ``(a/r)^2`` decay, the elliptical deformed hole) and it is exact on the hole boundary, which is where the finite element mesh has to work hardest.

**Plane stress** because the plate is thin: ``\sigma_z = \tau_{xz} = \tau_{yz} = 0``, and the thickness is free to change (``\varepsilon_z \neq 0``). The in-plane stresses of a traction problem are the same in plane stress and plane strain (Sadd §7.5), so the stress field below is also the stress field of a long cylinder with a hole; only the displacements would change (Table 7-1).
"""

# ╔═╡ e5000000-0000-0000-0000-000000000004
md"""
## 2. Governing equations in polar coordinates — plane stress

With no body forces, an Airy stress function ``\phi(r, \theta)`` generates stresses that satisfy the polar equilibrium equations identically (Sadd 7.6.7):

```math
\sigma_r = \frac{1}{r}\frac{\partial \phi}{\partial r} + \frac{1}{r^2}\frac{\partial^2 \phi}{\partial \theta^2}, \qquad
\sigma_\theta = \frac{\partial^2 \phi}{\partial r^2}, \qquad
\tau_{r\theta} = -\frac{\partial}{\partial r}\left(\frac{1}{r}\frac{\partial \phi}{\partial \theta}\right)
```

Compatibility is the biharmonic equation with the polar Laplacian (Sadd 7.6.8):

```math
\nabla^4 \phi = \left(\frac{\partial^2}{\partial r^2} + \frac{1}{r}\frac{\partial}{\partial r} + \frac{1}{r^2}\frac{\partial^2}{\partial \theta^2}\right)
\left(\frac{\partial^2 \phi}{\partial r^2} + \frac{1}{r}\frac{\partial \phi}{\partial r} + \frac{1}{r^2}\frac{\partial^2 \phi}{\partial \theta^2}\right) = 0
```

**Plane stress** Hooke's law in polar components (Sadd 7.2.3 and 7.6.2):

```math
\varepsilon_r = \frac{1}{E}\left(\sigma_r - \nu\sigma_\theta\right), \qquad
\varepsilon_\theta = \frac{1}{E}\left(\sigma_\theta - \nu\sigma_r\right), \qquad
\gamma_{r\theta} = \frac{2(1+\nu)}{E}\tau_{r\theta}, \qquad
\varepsilon_z = -\frac{\nu}{E}\left(\sigma_r + \sigma_\theta\right)
```

(compare Lecture 13, where plane strain had ``\varepsilon_z = 0`` and a non-zero ``\sigma_z`` instead). Strain–displacement relations in polar coordinates (Sadd 7.6.1):

```math
\varepsilon_r = \frac{\partial u_r}{\partial r}, \qquad
\varepsilon_\theta = \frac{1}{r}\left(u_r + \frac{\partial u_\theta}{\partial \theta}\right), \qquad
\gamma_{r\theta} = \frac{1}{r}\frac{\partial u_r}{\partial \theta} + \frac{\partial u_\theta}{\partial r} - \frac{u_\theta}{r}
```
"""

# ╔═╡ e5000000-0000-0000-0000-000000000005
md"""
## 3. Boundary conditions

**Hole** ``r = a``, for all ``\theta``: traction-free,

```math
\sigma_r(a, \theta) = 0, \qquad \tau_{r\theta}(a, \theta) = 0
```

**Far field** ``r \to \infty``: the stress must approach the uniform uniaxial state ``\sigma_x = T``, ``\sigma_y = \tau_{xy} = 0`` of Lecture 9. In polar components (stress transformation, Sadd 8.4.11):

```math
\sigma_r \to \frac{T}{2}\left(1 + \cos 2\theta\right), \qquad
\sigma_\theta \to \frac{T}{2}\left(1 - \cos 2\theta\right), \qquad
\tau_{r\theta} \to -\frac{T}{2}\sin 2\theta
```

We will not type these in: the uniform field has the Airy function ``\phi_0 = \tfrac{1}{2}T y^2`` (Lecture 9's ``A_{02} = T/2``), and with ``y = r\sin\theta`` that is ``\phi_0 = \tfrac{T}{4}r^2(1 - \cos 2\theta)``. Running it through the polar stress formulas gives the far-field stresses, and it also tells us which Michell terms to use: a ``r^2`` term and a ``r^2\cos 2\theta`` term.

**Symmetry.** The problem is symmetric about both the ``x`` and the ``y`` axes. The centre of the hole therefore does not move and nothing rotates: the rigid-body part of the displacement is zero.

### What is relaxed

* **Infinite plate.** Discussed in Section 1 and quantified in Section 9: the long edges are only ``2a`` from the hole centre.
* **Uniform end traction.** A testing machine grips the ends; the actual traction distribution there is not uniform. Saint-Venant (Sadd §5.6) confines that difference to within about one width ``W`` of each end, and the hole is ``3W`` away.
* **Plane stress.** Exact only for a vanishing thickness; for a plate of finite thickness ``h`` the hoop stress at the hole varies slightly through the thickness when ``h`` is comparable with ``a``. In the project the FEA model is two-dimensional plane stress, so this is not a source of difference.
"""

# ╔═╡ e5000000-0000-0000-0000-000000000006
md"""
## 4. Symbolic setup

Independent variables, problem parameters and the differential operators. The helpers are those of Lecture 13 (polar Laplacian, biharmonic, stresses, Laurent-polynomial integration in ``r``, exact Gauss–Jordan solver) plus two new ones for the ``\theta`` dependence.
"""

# ╔═╡ e5000000-0000-0000-0000-000000000007
begin
	@variables r θ
	@variables T a E ν             # far-field tension, hole radius, Young's modulus, Poisson's ratio
	Dr = Differential(r)
	Dθ = Differential(θ)
end

# ╔═╡ e5000000-0000-0000-0000-000000000008
"""
    laplacian_polar(f)

The two-dimensional Laplacian in polar coordinates, ∇² = ∂²/∂r² + (1/r) ∂/∂r + (1/r²) ∂²/∂θ² (Sadd 7.6.5).
"""
laplacian_polar(f) = expand_derivatives(Dr(Dr(f)) + Dr(f)/r + Dθ(Dθ(f))/r^2)

# ╔═╡ e5000000-0000-0000-0000-000000000009
"""
    biharmonic_polar(ϕ)

∇⁴ϕ = ∇²(∇²ϕ) in polar coordinates (Sadd 7.6.8).
"""
biharmonic_polar(ϕ) = simplify(laplacian_polar(laplacian_polar(ϕ)))

# ╔═╡ e5000000-0000-0000-0000-00000000000a
"""
    stresses_polar(ϕ) -> (σr, σθ, τrθ)

Stresses from the Airy stress function in polar coordinates (Sadd 7.6.7).
"""
stresses_polar(ϕ) = (simplify(expand_derivatives(Dr(ϕ)/r + Dθ(Dθ(ϕ))/r^2)),
                     simplify(expand_derivatives(Dr(Dr(ϕ)))),
                     simplify(expand_derivatives(-Dr(Dθ(ϕ)/r))))

# ╔═╡ e5000000-0000-0000-0000-00000000000b
md"""
Helpers for "for all ``r``" and "for all ``\theta``" conditions. Every stress in this problem has the form ``f_0(r) + f_2(r)\cos 2\theta + g_2(r)\sin 2\theta`` with ``f, g`` Laurent polynomials in ``r`` (powers ``r^{-4}, r^{-2}, r^0, r^2, \dots``):

* `laurent_coeffs(expr)` returns the coefficient of each power of ``r`` (the coefficients may still contain ``\theta``). A condition at ``r \to \infty`` is a condition on the ``r^0`` coefficient alone, because every negative power decays.
* `harmonic_coeffs(expr)` splits an expression into its ``1``, ``\cos 2\theta`` and ``\sin 2\theta`` parts. A condition that holds *for all* ``\theta`` is three conditions, one per part.
* `∫r` integrates a Laurent polynomial in ``r`` term by term (same as Lecture 13).
* `at_angle` evaluates at ``\theta = 0, \pi/2, \pi`` using exact values of the trigonometric functions, so no ``10^{-16}`` noise enters the exact row reduction.
"""

# ╔═╡ e5000000-0000-0000-0000-00000000000c
begin
	deg(mono, v) = mono isa Number ? 0 : Symbolics.degree(mono, v)

	function laurent_coeffs(expr; k = 4)
		e = expand(Symbolics.simplify_fractions(expand(expr) * r^k))
		d, rem = Symbolics.polynomial_coeffs(e, [r])
		isequal(Symbolics.value(rem), 0) || error("expression is not a Laurent polynomial in r")
		return Dict(deg(mono, r) - k => co for (mono, co) in d)
	end
	laurent_coeff(expr, m; k = 4) = get(laurent_coeffs(expr; k), m, 0)

	function harmonic_coeffs(expr)
		C, S = @variables C S                                                  # stand in for cos 2θ and sin 2θ
		e = expand(substitute(expand(expr), Dict(cos(2θ) => C, sin(2θ) => S)))
		d, rem = Symbolics.polynomial_coeffs(e, [C, S])
		isequal(Symbolics.value(rem), 0) || error("expression is not of the form f₀ + f₂ cos 2θ + g₂ sin 2θ")
		c0, c2, s2 = 0, 0, 0
		for (mono, co) in d
			dc, ds = deg(mono, C), deg(mono, S)
			if dc == 0 && ds == 0
				c0 += co
			elseif dc == 1 && ds == 0
				c2 += co
			elseif dc == 0 && ds == 1
				s2 += co
			else
				error("higher harmonics or products of cos 2θ and sin 2θ are not handled")
			end
		end
		return (c0 = c0, c2 = c2, s2 = s2)
	end

	function ∫r(expr; k = 4)
		ℓ = only(@variables ℓ)                                                # stands in for log(r) while we read coefficients
		e = expand(Symbolics.simplify_fractions(expand(expr) * r^k))
		e = expand(substitute(e, Dict(log(r) => ℓ)))
		d, rem = Symbolics.polynomial_coeffs(e, [r, ℓ])
		isequal(Symbolics.value(rem), 0) || error("expression is not a Laurent polynomial in r (times log r)")
		s = 0
		for (mono, co) in d
			m = deg(mono, r) - k
			n = deg(mono, ℓ)
			n ≤ 1 || error("powers of log r above the first are not handled")
			if n == 0
				s += m == -1 ? co*log(r) : co*r^(m + 1)/(m + 1)
			else
				s += m == -1 ? co*log(r)^2/2 : co*(r^(m + 1)/(m + 1)*log(r) - r^(m + 1)/(m + 1)^2)
			end
		end
		return expand(s)
	end
	∫r(expr, lo, hi; k = 4) = expand(substitute(∫r(expr; k), Dict(r => hi)) - substitute(∫r(expr; k), Dict(r => lo)))

	# exact evaluation at θ₀ given sin θ₀ and cos θ₀ (integers for θ₀ = 0, π/2, π)
	function at_angle(expr, θ₀, sinθ₀, cosθ₀)
		e = substitute(expr, Dict(sin(2θ) => 2sinθ₀*cosθ₀, cos(2θ) => cosθ₀^2 - sinθ₀^2))
		e = substitute(e, Dict(sin(θ) => sinθ₀, cos(θ) => cosθ₀))
		return substitute(e, Dict(θ => θ₀))
	end
end

# ╔═╡ e5000000-0000-0000-0000-00000000000d
"""
    solve_linear_bcs(eqs, unknowns) -> (solution::Dict, n_redundant)

Exact Gauss–Jordan elimination on the augmented matrix [A | b] of a linear system that may have
more equations than unknowns (same as Lectures 12–13). Numeric pivots are preferred, so we never divide by a
symbolic parameter unless we must. Errors if the system is inconsistent or leaves an unknown undetermined.
"""
function solve_linear_bcs(eqs, unknowns)
	isz(x)  = isequal(Symbolics.value(simplify(x)), 0)
	isnum(x) = Symbolics.value(x) isa Number
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

# ╔═╡ e5000000-0000-0000-0000-00000000000e
md"""
### The uniform field in polar coordinates

``\phi_0 = \tfrac{T}{4}r^2(1 - \cos 2\theta)`` is ``\tfrac12 T y^2``. Its stresses are the far-field values we have to match:
"""

# ╔═╡ e5000000-0000-0000-0000-00000000000f
begin
	ϕ0 = T*r^2*(1 - cos(2θ))/4
	σr_far, σθ_far, τrθ_far = stresses_polar(ϕ0)
	(σr_far = σr_far, σθ_far = σθ_far, τrθ_far = τrθ_far)
end

# ╔═╡ e5000000-0000-0000-0000-000000000010
md"""
## 5. Trial stress function: axisymmetric and ``\cos 2\theta`` terms of the Michell solution

The uniform field needs ``r^2`` and ``r^2\cos 2\theta``. The hole disturbs it, and we expect the disturbance to be *local*: to decay with distance from the hole. The Michell solution (Sadd 8.3.6) offers, for the same two angular dependences, exactly the terms with negative powers of ``r`` that we need. Following Sadd (8.4.12) we take every axisymmetric and every ``\cos 2\theta`` term:

```math
\phi = a_1\ln r + a_2 r^2 + a_3 r^2\ln r + \left(b_1 r^2 + b_2 r^4 + \frac{b_3}{r^2} + b_4\right)\cos 2\theta
```

(``a_0`` produces no stress and is dropped; ``\sin 2\theta`` terms are excluded by the symmetry about the ``x`` axis, and higher harmonics by the far field, which has none.) Sadd's ``a_{21}, a_{22}, a_{23}, a_{24}`` are our ``b_1, b_2, b_3, b_4``.
"""

# ╔═╡ e5000000-0000-0000-0000-000000000011
begin
	@variables a1 a2 a3 b1 b2 b3 b4
	ϕ_trial = a1*log(r) + a2*r^2 + a3*r^2*log(r) + (b1*r^2 + b2*r^4 + b3/r^2 + b4)*cos(2θ)
end

# ╔═╡ e5000000-0000-0000-0000-000000000012
md"**Biharmonic check.** Every Michell term is biharmonic by construction, so this should be zero with no condition on the coefficients:"

# ╔═╡ e5000000-0000-0000-0000-000000000013
∇⁴ϕ = biharmonic_polar(ϕ_trial)

# ╔═╡ e5000000-0000-0000-0000-000000000014
md"**Stresses from the trial function** (Sadd 8.4.13):"

# ╔═╡ e5000000-0000-0000-0000-000000000015
σr_trial, σθ_trial, τrθ_trial = stresses_polar(ϕ_trial)

# ╔═╡ e5000000-0000-0000-0000-000000000016
md"""
```math
\begin{aligned}
\sigma_r &= a_3(1 + 2\ln r) + 2a_2 + \frac{a_1}{r^2} - \left(2b_1 + \frac{6b_3}{r^4} + \frac{4b_4}{r^2}\right)\cos 2\theta \\
\sigma_\theta &= a_3(3 + 2\ln r) + 2a_2 - \frac{a_1}{r^2} + \left(2b_1 + 12 b_2 r^2 + \frac{6b_3}{r^4}\right)\cos 2\theta \\
\tau_{r\theta} &= \left(2b_1 + 6b_2 r^2 - \frac{6b_3}{r^4} - \frac{2b_4}{r^2}\right)\sin 2\theta
\end{aligned}
```

### Boundedness at infinity

Two terms grow without limit as ``r \to \infty``: ``a_3 r^2 \ln r`` gives stresses ``\propto \ln r`` and ``b_2 r^4`` gives stresses ``\propto r^2``. Neither can tend to the uniform field, so ``a_3 = b_2 = 0`` (Sadd: "for finite stresses at infinity"). The ``a_3`` term is doubly excluded: Lecture 13 showed that its tangential displacement grows linearly with ``\theta``, and a plate that *surrounds* the hole needs single-valued displacements.
"""

# ╔═╡ e5000000-0000-0000-0000-000000000017
begin
	bounded_map = Dict(a3 => 0, b2 => 0)
	σr_b  = substitute(σr_trial,  bounded_map)
	σθ_b  = substitute(σθ_trial,  bounded_map)
	τrθ_b = substitute(τrθ_trial, bounded_map)
	(σr_b, σθ_b, τrθ_b)
end

# ╔═╡ e5000000-0000-0000-0000-000000000018
md"""
## 6. What the far field and the hole say

**Far field.** The remaining terms are Laurent polynomials in ``r`` whose negative powers vanish at infinity, so the limit ``r \to \infty`` of each stress is its ``r^0`` coefficient. That coefficient must equal the uniform field *for every* ``\theta``, i.e. separately for the ``1``, ``\cos 2\theta`` and ``\sin 2\theta`` parts (``\sigma_\theta`` gives nothing new, but we include it and let the solver report it as redundant):
"""

# ╔═╡ e5000000-0000-0000-0000-000000000019
far_eqs = let
	eqs = Equation[]
	for (σ, σ∞) in ((σr_b, σr_far), (σθ_b, σθ_far), (τrθ_b, τrθ_far))
		h = harmonic_coeffs(laurent_coeff(σ - σ∞, 0))
		append!(eqs, [h.c0 ~ 0, h.c2 ~ 0, h.s2 ~ 0])
	end
	eqs
end

# ╔═╡ e5000000-0000-0000-0000-00000000001a
md"""
**Hole.** ``\sigma_r(a, \theta) = 0`` and ``\tau_{r\theta}(a, \theta) = 0`` for all ``\theta`` — again one equation per harmonic (Sadd 8.4.14):
"""

# ╔═╡ e5000000-0000-0000-0000-00000000001b
hole_eqs = let
	eqs = Equation[]
	for σ in (σr_b, τrθ_b)
		h = harmonic_coeffs(substitute(σ, Dict(r => a)))
		append!(eqs, [h.c0 ~ 0, h.c2 ~ 0, h.s2 ~ 0])
	end
	eqs
end

# ╔═╡ e5000000-0000-0000-0000-00000000001c
coeff_map, n_redundant = solve_linear_bcs(vcat(far_eqs, hole_eqs), [a1, a2, b1, b3, b4])

# ╔═╡ e5000000-0000-0000-0000-00000000001d
md"""
$(length(far_eqs) + length(hole_eqs)) equations (several of them ``0 = 0``), $(n_redundant) redundant, five unknowns, all determined:

```math
a_1 = -\frac{T a^2}{2}, \qquad a_2 = \frac{T}{4}, \qquad b_1 = -\frac{T}{4}, \qquad b_3 = -\frac{T a^4}{4}, \qquad b_4 = \frac{T a^2}{2}
```

(Sadd's values below (8.4.14).) The far field alone fixes ``a_2`` and ``b_1`` — they are just the uniform field ``\phi_0``. The hole then fixes the three decaying terms: ``a_1\ln r`` is the axisymmetric *pressurised-hole* field of Lecture 13 (here with the sign of a hole that is pulled open), and ``b_3/r^2`` and ``b_4`` are the two decaying ``\cos 2\theta`` terms, with stresses ``\propto 1/r^4`` and ``1/r^2``.

The stress function:
"""

# ╔═╡ e5000000-0000-0000-0000-00000000001e
ϕ = simplify(substitute(ϕ_trial, merge(bounded_map, coeff_map)))

# ╔═╡ e5000000-0000-0000-0000-00000000001f
md"""
```math
\phi = \frac{T}{4}\left[r^2 - 2a^2\ln r - \left(r^2 - 2a^2 + \frac{a^4}{r^2}\right)\cos 2\theta\right]
```

## 7. The stress field and the stress concentration factor
"""

# ╔═╡ e5000000-0000-0000-0000-000000000020
begin
	σr, σθ, τrθ = stresses_polar(ϕ)
	σz = 0                                   # plane stress
	(σr = σr, σθ = σθ, τrθ = τrθ)
end

# ╔═╡ e5000000-0000-0000-0000-000000000021
md"""
This is the Kirsch solution, Sadd (8.4.15):

```math
\begin{aligned}
\sigma_r &= \frac{T}{2}\left(1 - \frac{a^2}{r^2}\right) + \frac{T}{2}\left(1 - \frac{4a^2}{r^2} + \frac{3a^4}{r^4}\right)\cos 2\theta \\
\sigma_\theta &= \frac{T}{2}\left(1 + \frac{a^2}{r^2}\right) - \frac{T}{2}\left(1 + \frac{3a^4}{r^4}\right)\cos 2\theta \\
\tau_{r\theta} &= -\frac{T}{2}\left(1 + \frac{2a^2}{r^2} - \frac{3a^4}{r^4}\right)\sin 2\theta
\end{aligned}
```

Every term is the uniform field times a correction in powers of ``a/r``. Observations:

* **No ``E``, no ``\nu``, and no length other than ``a``.** The stresses depend on the material not at all and on the hole only through the ratio ``a/r``: a hole of any size in an infinite plate produces the same stress pattern, scaled. The same hole in a *finite* plate introduces a second length (``W``) and then the pattern does depend on ``a/W`` — Section 9.
* **The disturbance is local.** The hole's contribution decays like ``(a/r)^2`` at worst. Two hole radii away it is still a quarter of its peak; four radii away, a sixteenth.
* **On the hole** the only non-zero stress is the hoop stress (Sadd 8.4.16): ``\sigma_\theta(a, \theta) = T(1 - 2\cos 2\theta)``. It is ``3T`` at ``\theta = \pm 90°`` (the sides of the hole facing across the load), **zero** at ``\theta = \pm 30°, \pm 150°``, and ``-T`` at ``\theta = 0, 180°`` — a *compressive* stress on the hole surface, in a plate that is only pulled. The **stress concentration factor** is ``K_t = 3`` (Sadd 8.4.17).

Verify every condition we imposed (each residual should be zero):
"""

# ╔═╡ e5000000-0000-0000-0000-000000000022
bc_residuals = simplify.(expand.([
	biharmonic_polar(ϕ),
	substitute(σr,  Dict(r => a)),                         # hole free of normal traction, for all θ
	substitute(τrθ, Dict(r => a)),                         # hole free of shear traction, for all θ
	laurent_coeff(σr  - σr_far,  0),                       # far field, for all θ
	laurent_coeff(σθ  - σθ_far,  0),
	laurent_coeff(τrθ - τrθ_far, 0),
]))

# ╔═╡ e5000000-0000-0000-0000-000000000023
md"""
### The three lines that matter

The hoop stress **around the hole**, the normal stress **across the net section** (the ``y`` axis, ``\theta = \pi/2``, where ``\sigma_x = \sigma_\theta``) and the normal stress **along the load axis** (``\theta = 0``, where ``\sigma_y = \sigma_\theta``):
"""

# ╔═╡ e5000000-0000-0000-0000-000000000024
begin
	σθ_hole = simplify(substitute(σθ, Dict(r => a)))               # σθ(a, θ)
	σx_net  = simplify(at_angle(σθ, Num(π)/2, 1, 0))               # σx(0, y) with r = |y|
	σy_axis = simplify(at_angle(σθ, 0, 0, 1))                      # σy(x, 0) with r = |x|
	(σθ_hole = σθ_hole, σx_net = σx_net, σy_axis = σy_axis)
end

# ╔═╡ e5000000-0000-0000-0000-000000000025
md"""
```math
\sigma_\theta(a, \theta) = T\left(1 - 2\cos 2\theta\right), \qquad
\sigma_x(0, y) = T\left(1 + \frac{a^2}{2y^2} + \frac{3a^4}{2y^4}\right), \qquad
\sigma_y(x, 0) = \frac{T}{2}\left(\frac{a^2}{x^2} - \frac{3a^4}{x^4}\right)
```

The net-section formula is the one every stress-concentration handbook quotes. It is ``3T`` at ``y = a``, ``1.22\,T`` at ``y = 2a`` (where this strip's edge is) and ``1.07\,T`` at ``y = 3a``. The stress along the load axis is **compressive** right next to the hole, ``-T`` at ``x = a``, changes sign at ``x = \sqrt{3}\,a`` and never exceeds ``T/24`` in tension.

**Where does the load go?** In the infinite plate the net section is infinitely wide, but we can ask how much of the load passes through the part of it that *our strip actually has*, ``a \le |y| \le W/2``:
"""

# ╔═╡ e5000000-0000-0000-0000-000000000026
begin
	W = only(@variables W)
	F_net_strip = simplify(2 * ∫r(σx_net, a, W/2))     # ∫ σx(0, y) dy over both ligaments
	F_applied   = T*W                                                              # force per unit thickness at the ends
	(F_net_strip = F_net_strip, ratio = simplify(F_net_strip/F_applied), ratio_at_W_4a = simplify(substitute(F_net_strip/F_applied, Dict(W => 4a))))
end

# ╔═╡ e5000000-0000-0000-0000-000000000027
md"""
```math
\frac{1}{TW}\int_{\text{ligaments}}\sigma_x(0, y)\,dy = 1 - \frac{2a^2}{W^2}\left(1 + \frac{4a^2}{W^2}\right)
\qquad\Rightarrow\qquad \frac{27}{32} = 0.844 \;\text{ for } W = 4a
```

The infinite-plate stresses carry only 84 % of the applied load through the net section of our strip; the remaining 16 % is carried by material at ``|y| > W/2`` that the strip does not have. In the real strip that load has to squeeze through the ligaments too, so the true stresses on the net section are *higher* than Kirsch's — by roughly that factor on average, and by more at the hole edge. This is the finite-width effect in one line.

### Cartesian components

A finite element code reports ``\sigma_x, \sigma_y, \tau_{xy}``. Stress transformation from polar to Cartesian (Sadd 3.3, rotation by ``-\theta``):

```math
\sigma_x = \sigma_r\cos^2\theta + \sigma_\theta\sin^2\theta - 2\tau_{r\theta}\sin\theta\cos\theta, \quad
\sigma_y = \sigma_r\sin^2\theta + \sigma_\theta\cos^2\theta + 2\tau_{r\theta}\sin\theta\cos\theta, \quad
\tau_{xy} = (\sigma_r - \sigma_\theta)\sin\theta\cos\theta + \tau_{r\theta}(\cos^2\theta - \sin^2\theta)
```
"""

# ╔═╡ e5000000-0000-0000-0000-000000000028
begin
	σx  = simplify(σr*cos(θ)^2 + σθ*sin(θ)^2 - 2τrθ*sin(θ)*cos(θ))
	σy  = simplify(σr*sin(θ)^2 + σθ*cos(θ)^2 + 2τrθ*sin(θ)*cos(θ))
	τxy = simplify((σr - σθ)*sin(θ)*cos(θ) + τrθ*(cos(θ)^2 - sin(θ)^2))
	(σx = σx, σy = σy, τxy = τxy)
end

# ╔═╡ e5000000-0000-0000-0000-000000000029
md"""
Sanity check on the two symmetry lines, where the Cartesian and polar components coincide: on ``\theta = \pi/2`` we need ``\sigma_x = \sigma_\theta``, on ``\theta = 0`` we need ``\sigma_y = \sigma_\theta`` (both residuals zero):
"""

# ╔═╡ e5000000-0000-0000-0000-00000000002a
cartesian_check = simplify.([
	at_angle(σx, Num(π)/2, 1, 0) - σx_net,
	at_angle(σy, 0, 0, 1) - σy_axis,
])

# ╔═╡ e5000000-0000-0000-0000-00000000002b
md"""
## 8. Strains and displacements

Plane stress Hooke's law. The out-of-plane strain ``\varepsilon_z`` is not zero: the plate thins where ``\sigma_r + \sigma_\theta > 0`` and thickens where the sum is negative (next to the hole on the load axis).
"""

# ╔═╡ e5000000-0000-0000-0000-00000000002c
begin
	εr  = simplify((σr - ν*σθ)/E)
	εθ  = simplify((σθ - ν*σr)/E)
	γrθ = simplify(2(1 + ν)/E * τrθ)
	εz  = simplify(-ν/E * (σr + σθ))
	(εr = εr, εθ = εθ, γrθ = γrθ, εz = εz)
end

# ╔═╡ e5000000-0000-0000-0000-00000000002d
md"""
Integrate the strain–displacement relations, as in Lecture 13 but now with the ``\cos 2\theta`` dependence carried along.

**Step 1 — radial strain.** ``\partial u_r/\partial r = \varepsilon_r`` gives ``u_r = \int\varepsilon_r\,dr + f(\theta)``:
"""

# ╔═╡ e5000000-0000-0000-0000-00000000002e
ur_part = ∫r(εr)

# ╔═╡ e5000000-0000-0000-0000-00000000002f
md"""
**Step 2 — hoop strain.** ``\partial u_\theta/\partial\theta = r\varepsilon_\theta - \int\varepsilon_r\,dr - f(\theta)``. Split the known part ``h = r\varepsilon_\theta - \int\varepsilon_r\,dr`` into harmonics. Its ``\theta``-independent part must be zero — otherwise ``u_\theta`` would contain a term proportional to ``\theta`` and be multivalued around the hole (that was the ``a_3`` term, already gone) — and its ``\cos 2\theta`` part integrates to ``\tfrac12 h_2(r)\sin 2\theta``:
"""

# ╔═╡ e5000000-0000-0000-0000-000000000030
begin
	h_parts = harmonic_coeffs(r*εθ - ur_part)
	uθ_part = simplify(h_parts.c2/2) * sin(2θ)
	(h0 = simplify(h_parts.c0), h2 = simplify(h_parts.c2), uθ_part = uθ_part)
end

# ╔═╡ e5000000-0000-0000-0000-000000000031
md"""
**Step 3 — shear strain.** With ``u_r = u_r^{\text{part}} + f(\theta)`` and ``u_\theta = u_\theta^{\text{part}} - F(\theta) + g(r)``, ``F' = f``, the shear relation leaves (Lecture 13, Section 5) ``f = u_{x0}\cos\theta + u_{y0}\sin\theta`` and ``g = \omega_0 r``: a rigid translation and rotation, nothing else. The particular parts must already reproduce ``\gamma_{r\theta}`` on their own. Assemble the general field and verify all three strains (residuals zero):
"""

# ╔═╡ e5000000-0000-0000-0000-000000000032
begin
	@variables ux0 uy0 ω0
	ur_general = ur_part + ux0*cos(θ) + uy0*sin(θ)
	uθ_general = uθ_part - ux0*sin(θ) + uy0*cos(θ) + ω0*r
	(ur_general, uθ_general)
end

# ╔═╡ e5000000-0000-0000-0000-000000000033
strain_check_general = simplify.(expand.([
	expand_derivatives(Dr(ur_general)) - εr,
	expand_derivatives((ur_general + Dθ(uθ_general))/r) - εθ,
	expand_derivatives(Dθ(ur_general)/r + Dr(uθ_general) - uθ_general/r) - γrθ,
]))

# ╔═╡ e5000000-0000-0000-0000-000000000034
md"""
**Rigid-body motion.** Symmetry about both axes pins the hole centre and forbids rotation: ``u_{x0} = u_{y0} = \omega_0 = 0``. The displacement field is then
"""

# ╔═╡ e5000000-0000-0000-0000-000000000035
begin
	rigid_map = Dict(ux0 => 0, uy0 => 0, ω0 => 0)
	ur = simplify(substitute(ur_general, rigid_map))
	uθ = simplify(substitute(uθ_general, rigid_map))
	(ur = ur, uθ = uθ)
end

# ╔═╡ e5000000-0000-0000-0000-000000000036
md"""
```math
\begin{aligned}
u_r &= \frac{T}{2E}\left\{(1-\nu)\,r + (1+\nu)\frac{a^2}{r}
      + \left[(1+\nu)\,r + \frac{4a^2}{r} - (1+\nu)\frac{a^4}{r^3}\right]\cos 2\theta\right\} \\[4pt]
u_\theta &= -\frac{T}{2E}\left[(1+\nu)\,r + 2(1-\nu)\frac{a^2}{r} + (1+\nu)\frac{a^4}{r^3}\right]\sin 2\theta
\end{aligned}
```

The ``r``-proportional terms are the uniform field of Lecture 9 (``u = Tx/E``, ``v = -\nu Ty/E``) written in polar components; the ``1/r`` and ``1/r^3`` terms are the hole's contribution, which decays — but only like ``1/r``, one power slower than the stresses. The elasticity literature writes this with the shear modulus ``\mu = E/2(1+\nu)`` and Kolosov's constant ``\kappa = (3-\nu)/(1+\nu)`` (plane stress):

```math
u_r = \frac{Ta}{4\mu}\left[\frac{r}{a}\left(\frac{\kappa - 1}{2} + \cos 2\theta\right) + \frac{a}{r}\left(1 + (\kappa+1)\cos 2\theta\right) - \frac{a^3}{r^3}\cos 2\theta\right], \qquad
u_\theta = -\frac{Ta}{4\mu}\left[\frac{r}{a} + (\kappa - 1)\frac{a}{r} + \frac{a^3}{r^3}\right]\sin 2\theta
```

Check that ours agrees (both differences zero), and look at the hole itself:
"""

# ╔═╡ e5000000-0000-0000-0000-000000000037
begin
	μ = E/(2(1 + ν))
	κ = (3 - ν)/(1 + ν)
	ur_ref = T*a/(4μ) * ((r/a)*((κ - 1)/2 + cos(2θ)) + (a/r)*(1 + (κ + 1)*cos(2θ)) - (a^3/r^3)*cos(2θ))
	uθ_ref = -T*a/(4μ) * (r/a + (κ - 1)*(a/r) + a^3/r^3) * sin(2θ)
	kolosov_check = simplify.(expand.([ur - ur_ref, uθ - uθ_ref]))
end

# ╔═╡ e5000000-0000-0000-0000-000000000038
begin
	ur_hole   = simplify(substitute(ur, Dict(r => a)))                 # u_r(a, θ)
	ur_hole_x = simplify(at_angle(ur_hole, 0, 0, 1))                   # at θ = 0: along the load
	ur_hole_y = simplify(at_angle(ur_hole, Num(π)/2, 1, 0))            # at θ = π/2: across the load
	(ur_hole = ur_hole, ur_hole_x = ur_hole_x, ur_hole_y = ur_hole_y)
end

# ╔═╡ e5000000-0000-0000-0000-000000000039
md"""
```math
u_r(a, \theta) = \frac{Ta}{E}\left(1 + 2\cos 2\theta\right): \qquad u_r(a, 0) = \frac{3Ta}{E}, \qquad u_r(a, \pi/2) = -\frac{Ta}{E}
```

The hole becomes an **ellipse**: its diameter grows by ``6Ta/E`` along the load and shrinks by ``2Ta/E`` across it, and — remarkably — neither number depends on ``\nu`` (the ``\nu``-dependent parts of ``u_r`` cancel exactly at ``r = a``). For comparison, a line of length ``2a`` in the plate without a hole stretches by ``2Ta/E`` along the load and shrinks by ``2\nu Ta/E`` across it: the hole stretches three times as much as the material it replaced, the displacement version of ``K_t = 3``.

### Cartesian displacements

``u_x = u_r\cos\theta - u_\theta\sin\theta``, ``u_y = u_r\sin\theta + u_\theta\cos\theta`` — what the finite element code reports at its nodes:
"""

# ╔═╡ e5000000-0000-0000-0000-00000000003a
begin
	ux = simplify(ur*cos(θ) - uθ*sin(θ))
	uy = simplify(ur*sin(θ) + uθ*cos(θ))
	(ux = ux, uy = uy)
end

# ╔═╡ e5000000-0000-0000-0000-00000000003b
md"""
## 9. Numerical example — the plate in the figure

Steel, ``E = 200`` GPa, ``\nu = 0.3``; hole radius ``a = 12.5`` mm (diameter 25 mm), width ``W = 50`` mm, length ``L = 300`` mm. The figure gives no load, so we take ``T = 100`` MPa at the ends; every stress and displacement below is proportional to ``T``, so other loads are a rescaling. Units: N, mm, MPa.
"""

# ╔═╡ e5000000-0000-0000-0000-00000000003c
prm = (T = 100.0,        # MPa  far-field (end) tension
       a = 12.5,         # mm   hole radius (25 mm diameter)
       W = 50.0,         # mm   plate width
       L = 300.0,        # mm   plate length
       E = 200_000.0,    # MPa  steel
       ν = 0.3)

# ╔═╡ e5000000-0000-0000-0000-00000000003d
begin
	# build_function compiles each symbolic expression to a callable f(r, θ, T, a, E, ν)
	fargs = (r, θ, T, a, E, ν)
	compiled = Dict(name => build_function(expr, fargs...; expression = Val(false)) for (name, expr) in
		(:σr => σr, :σθ => σθ, :τrθ => τrθ, :σx => σx, :σy => σy, :τxy => τxy,
		 :ur => ur, :uθ => uθ, :ux => ux, :uy => uy, :εz => εz))

	# evaluate a field at polar (r, θ) or Cartesian (x, y); NaN inside the hole
	polar(name, ri, θi) = ri < prm.a*(1 - 1e-9) ? NaN : compiled[name](ri, θi, prm.T, prm.a, prm.E, prm.ν)   # tolerance keeps points on r = a
	cart(name, xi, yi)  = polar(name, hypot(xi, yi), atan(yi, xi))

	"""
	    kirsch(x, y) -> NamedTuple

	Every field of the infinite-plate (Kirsch) solution at the point (x, y), with the parameters `prm`.
	Stresses in MPa, displacements in mm, strain εz dimensionless. NaN inside the hole.
	"""
	kirsch(x, y) = (σx = cart(:σx, x, y), σy = cart(:σy, x, y), τxy = cart(:τxy, x, y),
	                σr = cart(:σr, x, y), σθ = cart(:σθ, x, y), τrθ = cart(:τrθ, x, y),
	                ux = cart(:ux, x, y), uy = cart(:uy, x, y), εz = cart(:εz, x, y))
	kirsch(0.0, prm.a)          # the top of the hole, where σx = σθ = 3T
end

# ╔═╡ e5000000-0000-0000-0000-00000000003e
md"""
### Key numbers

"Finite strip" in the last column is the accepted finite-width result discussed below; "—" means the infinite-plate value is all we have in closed form.
"""

# ╔═╡ e5000000-0000-0000-0000-00000000003f
begin
	# finite-width stress concentration factor for a circular hole in a strip under tension
	# (Pilkey & Pilkey, Peterson's Stress Concentration Factors, Chart 4.1 fit to Howland's 1930 solution):
	# K_tn is based on the net-section stress σ_net = T W/(W − d), K_tg on the gross stress T.
	Ktn_fit(d_W) = 3.000 - 3.140d_W + 3.667d_W^2 - 1.527d_W^3
	d_W   = 2prm.a/prm.W
	Ktn   = Ktn_fit(d_W)
	Ktg   = Ktn/(1 - d_W)
	σ_net = prm.T/(1 - d_W)
	(d_over_W = d_W, K_tn = Ktn, K_tg = Ktg, σ_net = σ_net, σmax_finite = Ktg*prm.T)
end

# ╔═╡ e5000000-0000-0000-0000-000000000040
key_numbers = let
	a_, T_, E_, W_ = prm.a, prm.T, prm.E, prm.W
	rows = [
		("σθ at the top/bottom of the hole (0, ±a)  [MPa]",    polar(:σθ, a_, π/2),      Ktg*T_),
		("σθ at the sides of the hole (±a, 0)  [MPa]",         polar(:σθ, a_, 0.0),      NaN),
		("σθ on the hole at θ = 30°  [MPa]",                    polar(:σθ, a_, π/6),      NaN),
		("σx at the strip edge on the net section (0, ±W/2)  [MPa]", cart(:σx, 0.0, W_/2), NaN),
		("mean σx on the net section  [MPa]",                   NaN,                      σ_net),
		("stress concentration factor on T",                    3.0,                      Ktg),
		("σy at the strip edge above the hole (0, ±W/2)  [MPa]", cart(:σy, 0.0, W_/2),    0.0),
		("σy at the sides of the hole (±a, 0)  [MPa]",          cart(:σy, a_, 0.0),       NaN),
		("hole diameter change along the load  [mm]",           2polar(:ur, a_, 0.0),     NaN),
		("hole diameter change across the load  [mm]",          2polar(:ur, a_, π/2),     NaN),
		("ux at the end (L/2, 0)  [mm]",                        cart(:ux, prm.L/2, 0.0),  NaN),
		("thickness strain εz at (0, ±a)  [10⁻⁶]",              1e6polar(:εz, a_, π/2),   NaN),
		("thickness strain εz at (±a, 0)  [10⁻⁶]",              1e6polar(:εz, a_, 0.0),   NaN),
	]
	fmt(v) = isnan(v) ? "—" : string(round(v, digits = 4) + 0.0)
	hdr = "| quantity | infinite plate (Kirsch) | finite strip |\n|---|---|---|\n"
	Markdown.parse(hdr * join(["| $(q) | $(fmt(p)) | $(fmt(f)) |" for (q, p, f) in rows], "\n"))
end

# ╔═╡ e5000000-0000-0000-0000-000000000041
md"""
With ``d/W = `` $(d_W) the net-section stress is already ``2T`` = $(σ_net) MPa, and Howland's finite-width solution puts the peak hoop stress at ``K_{tn} \sigma_{\text{net}}`` with ``K_{tn} \approx`` $(round(Ktn, digits = 2)), i.e. about $(round(Ktg, digits = 2)) ``T`` = $(round(Ktg*prm.T, digits = 0)) MPa against Kirsch's ``3T`` = $(3prm.T) MPa. **The infinite-plate solution underestimates the peak stress in this strip by about $(round(100*(Ktg/3 - 1), digits = 0)) %.** That is the single most important number for the project: a plane-stress finite element model with a fine mesh at the hole should land near $(round(Ktg*prm.T, digits = 0)) MPa, not 300, and a neural operator trained on finite-strip data must learn that ``K_t`` grows with ``a/W``. (Howland's series has ``d/W = 0.5`` as its practical limit; the fit above reproduces his tabulated values to about 1 %.)

### Stresses around the hole and along the symmetry lines
"""

# ╔═╡ e5000000-0000-0000-0000-000000000042
fig_lines = let
	fig = Figure(size = (1200, 400))
	T_, a_, W_ = prm.T, prm.a, prm.W

	ax1 = Axis(fig[1, 1], title = "Hoop stress on the hole boundary", xlabel = "θ [°]", ylabel = "σθ(a, θ) / T",
	           xticks = 0:45:360)
	θs = range(0, 2π, length = 361)
	lines!(ax1, rad2deg.(θs), polar.(:σθ, a_, θs) ./ T_, color = :firebrick, linewidth = 3, label = "Kirsch: 1 − 2cos 2θ")
	hlines!(ax1, [3.0], color = :firebrick, linestyle = :dot, linewidth = 1)
	hlines!(ax1, [Ktg], color = :black, linestyle = :dash, linewidth = 2, label = "finite strip peak, K_tg = $(round(Ktg, digits = 2))")
	hlines!(ax1, [0.0], color = (:gray, 0.5))
	axislegend(ax1, position = :rt)

	ax2 = Axis(fig[1, 2], title = "Across the net section, x = 0", xlabel = "y / a", ylabel = "σ / T")
	ys = range(a_, 4a_, length = 121)
	lines!(ax2, ys ./ a_, cart.(:σx, 0.0, ys) ./ T_, color = :firebrick, linewidth = 3, label = "σx (= σθ)")
	lines!(ax2, ys ./ a_, cart.(:σy, 0.0, ys) ./ T_, color = :steelblue, linewidth = 3, label = "σy (= σr)")
	vlines!(ax2, [W_/2a_], color = :black, linestyle = :dash, linewidth = 2, label = "strip edge y = W/2")
	hlines!(ax2, [1.0], color = (:gray, 0.5))
	axislegend(ax2, position = :rt)

	ax3 = Axis(fig[1, 3], title = "Along the load axis, y = 0", xlabel = "x / a", ylabel = "σ / T")
	xs = range(a_, 4a_, length = 121)
	lines!(ax3, xs ./ a_, cart.(:σx, xs, 0.0) ./ T_, color = :firebrick, linewidth = 3, label = "σx (= σr)")
	lines!(ax3, xs ./ a_, cart.(:σy, xs, 0.0) ./ T_, color = :steelblue, linewidth = 3, label = "σy (= σθ)")
	hlines!(ax3, [0.0, 1.0], color = (:gray, 0.5))
	axislegend(ax3, position = :rb)
	fig
end

# ╔═╡ e5000000-0000-0000-0000-000000000043
md"""
Left: the hoop stress around the hole, ``3T`` at the top and bottom, ``-T`` at the sides, zero at ``30°`` from the load axis; the dashed line is where the finite-strip peak sits. Middle: on the net section ``\sigma_x`` falls from ``3T`` at the hole to ``1.22\,T`` at the strip edge — it has *not* reached ``T`` when the strip runs out of material. ``\sigma_y`` is zero at the hole, peaks at ``0.375\,T`` at ``y = \sqrt2\,a`` and is still ``0.28\,T`` at the edge, where a free edge requires zero: this is the traction the real strip cannot supply. Right: along the load axis ``\sigma_x`` grows from zero at the hole back to ``T``, and ``\sigma_y`` is compressive (``-T`` at the hole) out to ``x = \sqrt3\,a``.

### Stress fields on the strip

``\sigma_x, \sigma_y, \tau_{xy}`` of the infinite-plate solution, painted on the region the strip occupies (the hole is blank). The upper row shows the whole 300 × 50 mm plate, the lower row a 100 × 50 mm window around the hole.
"""

# ╔═╡ e5000000-0000-0000-0000-000000000044
fig_fields = let
	T_, W_, L_ = prm.T, prm.W, prm.L
	fig = Figure(size = (1200, 760))
	fields = ((:σx, "σx / T", (-1.0, 3.0), :RdBu), (:σy, "σy / T", (-1.0, 1.0), :RdBu), (:τxy, "τxy / T", (-1.0, 1.0), :RdBu))
	for (j, (name, lab, lims, cmap)) in enumerate(fields)
		# whole strip
		xs = range(-L_/2, L_/2, length = 601); ys = range(-W_/2, W_/2, length = 101)
		Z = [cart(name, xi, yi)/T_ for xi in xs, yi in ys]
		ax = Axis(fig[j, 1], title = "$(lab), whole plate", xlabel = "x [mm]", ylabel = "y [mm]", aspect = DataAspect())
		heatmap!(ax, xs, ys, Z, colormap = cmap, colorrange = lims, nan_color = :white)
		lines!(ax, prm.a .* cos.(range(0, 2π, length = 181)), prm.a .* sin.(range(0, 2π, length = 181)), color = :black, linewidth = 1)
		# window around the hole
		xw = range(-2W_, 2W_, length = 401); yw = range(-W_/2, W_/2, length = 101)
		Zw = [cart(name, xi, yi)/T_ for xi in xw, yi in yw]
		axw = Axis(fig[j, 2], title = "$(lab), near the hole", xlabel = "x [mm]", ylabel = "y [mm]", aspect = DataAspect())
		hm = heatmap!(axw, xw, yw, Zw, colormap = cmap, colorrange = lims, nan_color = :white)
		lines!(axw, prm.a .* cos.(range(0, 2π, length = 181)), prm.a .* sin.(range(0, 2π, length = 181)), color = :black, linewidth = 1)
		Colorbar(fig[j, 3], hm, label = lab)
	end
	colsize!(fig.layout, 1, Relative(0.62))
	fig
end

# ╔═╡ e5000000-0000-0000-0000-000000000045
md"""
``\sigma_x`` is ``T`` almost everywhere — the hole's influence is confined to about two radii — with the ``3T`` hot spots at the top and bottom of the hole and the ``\sigma_x = 0`` shadows at its sides. ``\sigma_y`` shows the four-lobed ``\cos 2\theta`` pattern: compressive at the sides of the hole, tensile above and below it — and visibly non-zero at the strip edges ``y = \pm 25`` mm. ``\tau_{xy}`` has the antisymmetric ``\sin 2\theta`` pattern and is also non-zero on the edges. These edge values are the price of the infinite-plate idealisation; they are plotted next.

### What the infinite-plate solution leaves on the strip boundary

On the free edges ``y = \pm W/2`` the real plate has ``\sigma_y = \tau_{xy} = 0``; on the ends ``x = \pm L/2`` it has ``\sigma_x = T``, ``\tau_{xy} = 0``. The Kirsch field's values there are the **residual tractions** that a finite element model removes:
"""

# ╔═╡ e5000000-0000-0000-0000-000000000046
fig_edges = let
	T_, W_, L_ = prm.T, prm.W, prm.L
	fig = Figure(size = (1100, 380))
	ax1 = Axis(fig[1, 1], title = "Residual tractions on the long edge y = W/2", xlabel = "x [mm]", ylabel = "traction / T")
	xs = range(-L_/2, L_/2, length = 601)
	lines!(ax1, xs, cart.(:σy, xs, W_/2) ./ T_, color = :steelblue, linewidth = 3, label = "σy / T  (should be 0)")
	lines!(ax1, xs, cart.(:τxy, xs, W_/2) ./ T_, color = :darkorange, linewidth = 3, label = "τxy / T  (should be 0)")
	hlines!(ax1, [0.0], color = (:gray, 0.5))
	axislegend(ax1, position = :rt)

	ax2 = Axis(fig[1, 2], title = "Residual tractions on the end x = L/2", xlabel = "y [mm]", ylabel = "traction / T")
	ys = range(-W_/2, W_/2, length = 201)
	lines!(ax2, ys, (cart.(:σx, L_/2, ys) .- T_) ./ T_, color = :firebrick, linewidth = 3, label = "(σx − T) / T  (should be 0)")
	lines!(ax2, ys, cart.(:τxy, L_/2, ys) ./ T_, color = :darkorange, linewidth = 3, label = "τxy / T  (should be 0)")
	hlines!(ax2, [0.0], color = (:gray, 0.5))
	axislegend(ax2, position = :rt)
	fig
end

# ╔═╡ e5000000-0000-0000-0000-000000000047
md"""
At the **ends** the residuals are below $(round(100*maximum(abs.((cart.(:σx, prm.L/2, range(-prm.W/2, prm.W/2, length = 201)) .- prm.T) ./ prm.T)), digits = 2)) % of ``T``: the ends really are at infinity, and the 300 mm length could be halved without changing anything near the hole. At the **long edges** the residual normal traction reaches $(round(100*maximum(abs.(cart.(:σy, range(-prm.L/2, prm.L/2, length = 601), prm.W/2) ./ prm.T)), digits = 0)) % of ``T`` right above the hole and the shear reaches $(round(100*maximum(abs.(cart.(:τxy, range(-prm.L/2, prm.L/2, length = 601), prm.W/2) ./ prm.T)), digits = 0)) % of ``T`` at ``x \approx \pm a``, over a length of a few hole radii. Howland's correction is the elasticity solution of the strip loaded by the *negatives* of these tractions on its edges (and nothing on the hole); its hoop stress at the hole is what lifts ``K_t`` from 3 to about $(round(Ktg, digits = 1)). In the project, the FEA model does that superposition implicitly.

### Deformed shape

Displacements magnified. The hole turns into an ellipse elongated along the load; the plate narrows (Poisson) and the narrowing is strongest on the net section.
"""

# ╔═╡ e5000000-0000-0000-0000-000000000048
fig_deformed = let
	a_, W_, L_ = prm.a, prm.W, prm.L
	scale = 200
	fig = Figure(size = (1200, 560))

	# helper: a polyline (xs, ys) mapped to its deformed position, with points inside the hole removed
	function deformed(xs, ys, s)
		px = Float64[]; py = Float64[]
		for (xi, yi) in zip(xs, ys)
			if hypot(xi, yi) < a_*(1 - 1e-9)
				push!(px, NaN); push!(py, NaN)
			else
				push!(px, xi + s*cart(:ux, xi, yi)); push!(py, yi + s*cart(:uy, xi, yi))
			end
		end
		return px, py
	end
	function draw!(ax, xlim, dx, dy, s, col, lw)
		for x0 in range(-xlim, xlim, step = dx)
			ys = range(-W_/2, W_/2, length = 101)
			lines!(ax, deformed(fill(x0, length(ys)), ys, s)..., color = col, linewidth = lw)
		end
		for y0 in range(-W_/2, W_/2, step = dy)
			xs = range(-xlim, xlim, length = 401)
			lines!(ax, deformed(xs, fill(y0, length(xs)), s)..., color = col, linewidth = lw)
		end
		θs = range(0, 2π, length = 181)
		lines!(ax, deformed(a_ .* cos.(θs), a_ .* sin.(θs), s)..., color = col, linewidth = 2lw)
	end

	ax1 = Axis(fig[1, 1:2], title = "Whole plate, displacements × $scale", xlabel = "x [mm]", ylabel = "y [mm]", aspect = DataAspect())
	draw!(ax1, L_/2, 12.5, 12.5, 0, (:gray, 0.6), 1)
	draw!(ax1, L_/2, 12.5, 12.5, scale, :firebrick, 1)

	ax2 = Axis(fig[2, 1], title = "Near the hole, displacements × $scale", xlabel = "x [mm]", ylabel = "y [mm]", aspect = DataAspect())
	draw!(ax2, 2a_, 6.25, 6.25, 0, (:gray, 0.6), 1)
	draw!(ax2, 2a_, 6.25, 6.25, scale, :firebrick, 1.5)
	xlims!(ax2, -2.2a_, 2.2a_)

	ax3 = Axis(fig[2, 2], title = "Hole boundary, displacements × $scale", xlabel = "x [mm]", ylabel = "y [mm]", aspect = DataAspect())
	θs = range(0, 2π, length = 181)
	lines!(ax3, a_ .* cos.(θs), a_ .* sin.(θs), color = (:gray, 0.6), linewidth = 2, label = "undeformed circle")
	lines!(ax3, deformed(a_ .* cos.(θs), a_ .* sin.(θs), scale)..., color = :firebrick, linewidth = 3, label = "deformed (ellipse)")
	Legend(fig[2, 3], ax3)
	fig
end

# ╔═╡ e5000000-0000-0000-0000-000000000049
md"""
The hole's diameter grows by $(round(2polar(:ur, prm.a, 0.0)*1000, digits = 1)) μm along the load and shrinks by $(round(-2polar(:ur, prm.a, π/2)*1000, digits = 1)) μm across it, against a $(round(2prm.T*prm.a/prm.E*1000, digits = 1)) μm stretch of the same 25 mm of solid plate. The far ends move $(round(cart(:ux, prm.L/2, 0.0)*1000, digits = 0)) μm each — almost exactly the ``TL/2E`` of the plate without a hole, because the hole's extra compliance is local.

## 10. Using this solution in the project

`kirsch(x, y)` evaluates every field at a point, so the solution can be sampled on any grid — for example the nodes of the finite element mesh, or the pixel grid of a neural-operator input. A grid of the whole strip:
"""

# ╔═╡ e5000000-0000-0000-0000-00000000004a
grid_sample = let
	xs = range(-prm.L/2, prm.L/2, length = 121)
	ys = range(-prm.W/2, prm.W/2, length = 21)
	S = [kirsch(xi, yi) for xi in xs, yi in ys]              # 121 × 21 matrix of NamedTuples
	(size = size(S), max_σx_over_T = maximum(filter(!isnan, [s.σx for s in S]))/prm.T,
	 min_σy_over_T = minimum(filter(!isnan, [s.σy for s in S]))/prm.T)
end

# ╔═╡ e5000000-0000-0000-0000-00000000004b
md"""
(This grid happens to contain the points ``(0, \pm a)``, so it catches the ``3T`` peak. A grid, or a finite element mesh, with no node exactly there reports a lower peak — compare peak values on the hole boundary itself, not on a grid.)

Checks worth making against the finite element model of the 300 × 50 × ⌀25 plate, in order of usefulness:

1. **Hoop stress on the hole boundary**, ``\sigma_\theta(a, \theta)``. Kirsch gives ``T(1 - 2\cos 2\theta)``; the strip should show the same shape but a higher peak, ``\approx`` $(round(Ktg, digits = 1)) ``T`` instead of ``3T``, and a slightly different compressive value at ``\theta = 0``. The hole must be traction-free in both: ``\sigma_r = \tau_{r\theta} = 0``.
2. **``\sigma_x`` across the net section** ``x = 0``. The FEA curve must integrate to the applied force ``TW`` over the two ligaments; the Kirsch curve integrates to only ``0.84\,TW`` (Section 7), which is why the FEA curve lies above it everywhere.
3. **Compressive ``\sigma_y`` at the sides of the hole**, ``-T`` at ``(\pm a, 0)``. A mesh too coarse at the hole misses this as badly as it misses the ``3T``.
4. **Deformed hole**: elongation ``6Ta/E`` and contraction ``2Ta/E``, independent of ``\nu``. Finite width increases both a little.
5. **Far from the hole** (``|x| > 3a``) the FEA and Kirsch fields should agree to better than 1 % in ``\sigma_x``; if they do not, the end boundary conditions of the model are wrong, not the mesh.
6. **Edge tractions**: the FEA has ``\sigma_y = \tau_{xy} = 0`` on ``y = \pm W/2`` by construction; the difference FEA − Kirsch near the edge is Howland's correction field, and its decay away from the hole shows how far the finite-width effect reaches.

For the neural operator the same solution is a cheap source of *exact* data on the infinite-plate family (any ``a``, any ``T``, any ``E, \nu``) and a stress test: an operator that has learned the strip must reproduce the ``(a/r)^2`` decay and the ``\cos 2\theta`` structure that this notebook derives, and must deviate from them in the way Howland's correction does as ``a/W`` grows.
"""

# ╔═╡ e5000000-0000-0000-0000-00000000004c
md"""
## 11. Summary

| quantity | symbolic (infinite plate) | numeric, ``T`` = $(prm.T) MPa, ``a`` = $(prm.a) mm |
|---|---|---|
| ``\phi`` | ``\dfrac{T}{4}\left[r^2 - 2a^2\ln r - \left(r^2 - 2a^2 + \dfrac{a^4}{r^2}\right)\cos 2\theta\right]`` | — |
| ``\sigma_r`` | ``\dfrac{T}{2}\left(1 - \dfrac{a^2}{r^2}\right) + \dfrac{T}{2}\left(1 - \dfrac{4a^2}{r^2} + \dfrac{3a^4}{r^4}\right)\cos 2\theta`` | 0 on the hole |
| ``\sigma_\theta`` | ``\dfrac{T}{2}\left(1 + \dfrac{a^2}{r^2}\right) - \dfrac{T}{2}\left(1 + \dfrac{3a^4}{r^4}\right)\cos 2\theta`` | $(round(polar(:σθ, prm.a, π/2), digits = 1)) MPa at (0, ±a), $(round(polar(:σθ, prm.a, 0.0), digits = 1)) MPa at (±a, 0) |
| ``\tau_{r\theta}`` | ``-\dfrac{T}{2}\left(1 + \dfrac{2a^2}{r^2} - \dfrac{3a^4}{r^4}\right)\sin 2\theta`` | 0 on the hole |
| hoop stress on the hole | ``T(1 - 2\cos 2\theta)``, ``K_t = 3`` | peak $(3prm.T) MPa |
| net section ``\sigma_x(0, y)`` | ``T\left(1 + \dfrac{a^2}{2y^2} + \dfrac{3a^4}{2y^4}\right)`` | $(round(cart(:σx, 0.0, prm.W/2), digits = 1)) MPa at the strip edge |
| ``u_r`` on the hole | ``\dfrac{Ta}{E}(1 + 2\cos 2\theta)`` | diameter +$(round(2polar(:ur, prm.a, 0.0)*1000, digits = 1)) μm along, −$(round(-2polar(:ur, prm.a, π/2)*1000, digits = 1)) μm across the load |
| finite strip, ``d/W`` = $(d_W) | ``K_{tn} \approx 3 - 3.14\,\tfrac{d}{W} + 3.67\,\tfrac{d^2}{W^2} - 1.53\,\tfrac{d^3}{W^3}`` (Howland / Peterson) | ``K_{tg}`` = $(round(Ktg, digits = 2)), peak ≈ $(round(Ktg*prm.T, digits = 0)) MPa |

- The **polar Airy function** solves the hole exactly: of the Michell terms ``a_1\ln r + a_2 r^2 + a_3 r^2\ln r + (b_1 r^2 + b_2 r^4 + b_3 r^{-2} + b_4)\cos 2\theta``, boundedness at infinity kills ``a_3, b_2``, the far field fixes ``a_2, b_1`` (the uniform field ``\tfrac12 Ty^2``) and the traction-free hole fixes the three decaying terms ``a_1, b_3, b_4``.
- The stresses are **independent of ``E, \nu`` and of the hole size**; the hoop stress on the hole is ``T(1 - 2\cos 2\theta)``: ``3T`` at the top and bottom, ``-T`` (compression) at the sides. The disturbance decays like ``(a/r)^2``, the displacement perturbation like ``a/r``.
- The hole deforms into an **ellipse** with axis changes ``6Ta/E`` and ``-2Ta/E``, independent of ``\nu``.
- For the plate in the figure the long edges are only **two radii** from the hole. The infinite-plate field leaves tractions of up to ~20 % of ``T`` on them and carries only 84 % of the load through the net section; the finite-width solution (Howland) raises the peak stress from ``3T`` to about $(round(Ktg, digits = 1))``T``. The ends, twelve radii away, are at infinity to within 1 %.
- Everything a finite element model or a neural operator produces for this plate can be checked against `kirsch(x, y)` far from the edges, and must *differ* from it near the hole by the finite-width correction.
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
# ╠═e5000000-0000-0000-0000-000000000001
# ╟─e5000000-0000-0000-0000-000000000002
# ╟─e5000000-0000-0000-0000-000000000003
# ╟─e5000000-0000-0000-0000-000000000004
# ╟─e5000000-0000-0000-0000-000000000005
# ╟─e5000000-0000-0000-0000-000000000006
# ╠═e5000000-0000-0000-0000-000000000007
# ╠═e5000000-0000-0000-0000-000000000008
# ╠═e5000000-0000-0000-0000-000000000009
# ╠═e5000000-0000-0000-0000-00000000000a
# ╟─e5000000-0000-0000-0000-00000000000b
# ╠═e5000000-0000-0000-0000-00000000000c
# ╠═e5000000-0000-0000-0000-00000000000d
# ╟─e5000000-0000-0000-0000-00000000000e
# ╠═e5000000-0000-0000-0000-00000000000f
# ╟─e5000000-0000-0000-0000-000000000010
# ╠═e5000000-0000-0000-0000-000000000011
# ╟─e5000000-0000-0000-0000-000000000012
# ╠═e5000000-0000-0000-0000-000000000013
# ╟─e5000000-0000-0000-0000-000000000014
# ╠═e5000000-0000-0000-0000-000000000015
# ╟─e5000000-0000-0000-0000-000000000016
# ╠═e5000000-0000-0000-0000-000000000017
# ╟─e5000000-0000-0000-0000-000000000018
# ╠═e5000000-0000-0000-0000-000000000019
# ╟─e5000000-0000-0000-0000-00000000001a
# ╠═e5000000-0000-0000-0000-00000000001b
# ╠═e5000000-0000-0000-0000-00000000001c
# ╟─e5000000-0000-0000-0000-00000000001d
# ╠═e5000000-0000-0000-0000-00000000001e
# ╟─e5000000-0000-0000-0000-00000000001f
# ╠═e5000000-0000-0000-0000-000000000020
# ╟─e5000000-0000-0000-0000-000000000021
# ╠═e5000000-0000-0000-0000-000000000022
# ╟─e5000000-0000-0000-0000-000000000023
# ╠═e5000000-0000-0000-0000-000000000024
# ╟─e5000000-0000-0000-0000-000000000025
# ╠═e5000000-0000-0000-0000-000000000026
# ╟─e5000000-0000-0000-0000-000000000027
# ╠═e5000000-0000-0000-0000-000000000028
# ╟─e5000000-0000-0000-0000-000000000029
# ╠═e5000000-0000-0000-0000-00000000002a
# ╟─e5000000-0000-0000-0000-00000000002b
# ╠═e5000000-0000-0000-0000-00000000002c
# ╟─e5000000-0000-0000-0000-00000000002d
# ╠═e5000000-0000-0000-0000-00000000002e
# ╟─e5000000-0000-0000-0000-00000000002f
# ╠═e5000000-0000-0000-0000-000000000030
# ╟─e5000000-0000-0000-0000-000000000031
# ╠═e5000000-0000-0000-0000-000000000032
# ╠═e5000000-0000-0000-0000-000000000033
# ╟─e5000000-0000-0000-0000-000000000034
# ╠═e5000000-0000-0000-0000-000000000035
# ╟─e5000000-0000-0000-0000-000000000036
# ╠═e5000000-0000-0000-0000-000000000037
# ╠═e5000000-0000-0000-0000-000000000038
# ╟─e5000000-0000-0000-0000-000000000039
# ╠═e5000000-0000-0000-0000-00000000003a
# ╟─e5000000-0000-0000-0000-00000000003b
# ╠═e5000000-0000-0000-0000-00000000003c
# ╠═e5000000-0000-0000-0000-00000000003d
# ╟─e5000000-0000-0000-0000-00000000003e
# ╠═e5000000-0000-0000-0000-00000000003f
# ╠═e5000000-0000-0000-0000-000000000040
# ╟─e5000000-0000-0000-0000-000000000041
# ╠═e5000000-0000-0000-0000-000000000042
# ╟─e5000000-0000-0000-0000-000000000043
# ╠═e5000000-0000-0000-0000-000000000044
# ╟─e5000000-0000-0000-0000-000000000045
# ╠═e5000000-0000-0000-0000-000000000046
# ╟─e5000000-0000-0000-0000-000000000047
# ╠═e5000000-0000-0000-0000-000000000048
# ╟─e5000000-0000-0000-0000-000000000049
# ╠═e5000000-0000-0000-0000-00000000004a
# ╟─e5000000-0000-0000-0000-00000000004b
# ╟─e5000000-0000-0000-0000-00000000004c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
