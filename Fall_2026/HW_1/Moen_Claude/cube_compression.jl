# cube_compression.jl
#
# Deformation of a linear-elastic 3D cube under a uniform compressive stress on
# its top surface, with the bottom surface restrained.  Solved with Ferrite.jl.
#
#          z (w)
#          ^        σ (uniform, pressing down on the top face)
#          |   | | | | | |
#          |   v v v v v v
#          +-----------------+
#         /                 /|
#        +-----------------+ |  H = 28 mm
#        |                 | |
#        |                 | +   ---> x (u)
#        |                 |/
#        +-----------------+
#        //////////////////   bottom restrained
#      y (v)
#
# Units: mm, N, MPa (N/mm^2).
#
# Two bottom boundary conditions are studied:
#   :clamped     u = v = w = 0 on the bottom  (bottom glued to a rigid platen)
#   :frictionless w = 0 on the bottom, lateral sliding allowed (rigid-body modes
#                 removed by u = 0 on x = H/2 and v = 0 on y = H/2).  The exact answer is homogeneous uniaxial
#                 stress, so this case verifies the code.
#   :rough_platens bottom clamped AND the top pressed by a rigid, rough platen:
#                 u = v = 0 on the top face and w = -Δ (uniform).  Δ is scaled
#                 after the solve so the platen force equals σ0 H² (linear problem).
#
# Run:  julia --project=. cube_compression.jl

using Ferrite, SparseArrays, Printf, LinearAlgebra, Logging

# ---------------------------------------------------------------------------
# Problem data
# ---------------------------------------------------------------------------
const H  = 28.0        # cube edge length [mm]
# Material: elastic thermoplastic rubber (TPR / TPE, e.g. SEBS or Santoprene-type,
# Shore A ≈ 60-70).  Rubbers are nearly incompressible; ν is kept at 0.49 rather
# than 0.5 because ν -> 0.5 makes the linear-hex displacement formulation lock
# (λ/μ = 2ν/(1-2ν) = 49 here; it would be infinite at ν = 0.5).
const E  = 4.0         # Young's modulus [MPa]
const ν  = 0.49        # Poisson's ratio
# Load: a 5 oz can of tuna resting on the top face.  The can (≈85 mm diameter)
# is wider than the cube, so its whole weight is carried by the 28 mm x 28 mm face.
const m_can = 5.0 * 0.028349523    # 5 oz -> kg
const g     = 9.81                 # m/s^2
const W_can = m_can * g            # weight [N]
const σ0    = W_can / H^2          # uniform pressure on the top face [MPa]  (≈ 1.77e-3 MPa = 1.77 kPa)
const nel   = 14       # elements per edge (2 mm hexahedra)
const order = 2        # displacement interpolation order.  Use 2: with ν = 0.49 the
                       # linear hexahedron shows spurious stress oscillations at the
                       # clamped faces and converges slowly (volumetric locking).

# ---------------------------------------------------------------------------
# Elasticity tensor  σ = C : ε   (isotropic, linear elastic)
# ---------------------------------------------------------------------------
function elasticity_tensor(E, ν)
    λ = E * ν / ((1 + ν) * (1 - 2ν))
    μ = E / (2(1 + ν))
    δ(i, j) = i == j ? 1.0 : 0.0
    return SymmetricTensor{4,3}((i, j, k, l) ->
        λ * δ(i, j) * δ(k, l) + μ * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k)))
end

# ---------------------------------------------------------------------------
# Element stiffness  ke_ij = ∫ ∇ˢNi : C : ∇ˢNj dΩ
# ---------------------------------------------------------------------------
function assemble_cell!(ke, cv, C)
    n = getnbasefunctions(cv)
    for q in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q)
        for i in 1:n
            ∇Ni = shape_symmetric_gradient(cv, q, i)
            ∇Ni_C = ∇Ni ⊡ C
            for j in 1:n
                ∇Nj = shape_symmetric_gradient(cv, q, j)
                ke[i, j] += (∇Ni_C ⊡ ∇Nj) * dΩ
            end
        end
    end
    return ke
end

function assemble_stiffness!(K, dh, cv, C)
    n = getnbasefunctions(cv)
    ke = zeros(n, n)
    assembler = start_assemble(K)
    for cell in CellIterator(dh)
        reinit!(cv, cell)
        fill!(ke, 0.0)
        assemble_cell!(ke, cv, C)
        assemble!(assembler, celldofs(cell), ke)
    end
    return K
end

# ---------------------------------------------------------------------------
# Surface traction  f_i = ∫_Γ Ni ⋅ t dΓ   on a facet set
# ---------------------------------------------------------------------------
function assemble_traction!(f, dh, fv, facetset, t::Vec{3})
    n = getnbasefunctions(fv)
    fe = zeros(n)
    for facet in FacetIterator(dh, facetset)
        reinit!(fv, facet)
        fill!(fe, 0.0)
        for q in 1:getnquadpoints(fv)
            dΓ = getdetJdV(fv, q)
            for i in 1:n
                fe[i] += (shape_value(fv, q, i) ⋅ t) * dΓ
            end
        end
        assemble!(f, celldofs(facet), fe)
    end
    return f
end

# ---------------------------------------------------------------------------
# Strain and stress recovery at quadrature points, then L2-projection to nodes
# ---------------------------------------------------------------------------
function compute_qp_strain_stress(dh, cv, u, C)
    T = SymmetricTensor{2,3,Float64,6}
    ε_qp = Vector{Vector{T}}(undef, getncells(dh.grid))
    σ_qp = Vector{Vector{T}}(undef, getncells(dh.grid))
    for cell in CellIterator(dh)
        reinit!(cv, cell)
        ue = u[celldofs(cell)]
        εc = T[]; σc = T[]
        for q in 1:getnquadpoints(cv)
            ε = function_symmetric_gradient(cv, q, ue)
            push!(εc, ε)
            push!(σc, C ⊡ ε)
        end
        ε_qp[cellid(cell)] = εc
        σ_qp[cellid(cell)] = σc
    end
    return ε_qp, σ_qp
end

# ---------------------------------------------------------------------------
# Solve one case
# ---------------------------------------------------------------------------
function solve_cube(bottom_bc::Symbol; nel = nel, order = order, verbose = true)
    grid = generate_grid(Hexahedron, (nel, nel, nel), Vec(0.0, 0.0, 0.0), Vec(H, H, H))

    ip  = Lagrange{RefHexahedron,order}()^3      # vector-valued Lagrange interpolation
    qr  = QuadratureRule{RefHexahedron}(order + 1)   # full Gauss integration
    fqr = FacetQuadratureRule{RefHexahedron}(order + 1)
    cv  = CellValues(qr, ip)
    fv  = FacetValues(fqr, ip)

    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    # --- boundary conditions -----------------------------------------------
    ch = ConstraintHandler(dh)
    bottom = getfacetset(grid, "bottom")
    top    = getfacetset(grid, "top")
    Δtrial = σ0 * H / E               # trial platen displacement, rescaled after solve
    if bottom_bc == :clamped
        add!(ch, Dirichlet(:u, bottom, (x, t) -> Vec(0.0, 0.0, 0.0), [1, 2, 3]))
    elseif bottom_bc == :rough_platens
        add!(ch, Dirichlet(:u, bottom, (x, t) -> Vec(0.0, 0.0, 0.0), [1, 2, 3]))
        add!(ch, Dirichlet(:u, top,    (x, t) -> Vec(0.0, 0.0, -Δtrial), [1, 2, 3]))
    elseif bottom_bc == :frictionless
        add!(ch, Dirichlet(:u, bottom, (x, t) -> 0.0, [3]))     # w = 0 only
        # Remove the in-plane rigid-body modes symmetrically: u = 0 on the plane
        # x = H/2 and v = 0 on the plane y = H/2 (needs even nel so the planes carry
        # nodes).  This keeps the frame centred like the clamped cases and does not
        # disturb the homogeneous strain field.  Only vertex dofs are pinned, which
        # is all that is needed; Ferrite warns about nodeset constraints with
        # quadratic interpolation, so that warning is silenced here.
        iseven(nel) || error("frictionless case needs an even nel")
        addnodeset!(grid, "xmid", x -> abs(x[1] - H / 2) < 1e-8)
        addnodeset!(grid, "ymid", x -> abs(x[2] - H / 2) < 1e-8)
        Logging.with_logger(Logging.SimpleLogger(stderr, Logging.Error)) do
            add!(ch, Dirichlet(:u, getnodeset(grid, "xmid"), (x, t) -> 0.0, [1]))
            add!(ch, Dirichlet(:u, getnodeset(grid, "ymid"), (x, t) -> 0.0, [2]))
        end
    else
        error("unknown bottom_bc = $bottom_bc")
    end
    close!(ch)

    # --- assemble ------------------------------------------------------------
    C = elasticity_tensor(E, ν)
    K = allocate_matrix(dh)
    f = zeros(ndofs(dh))
    assemble_stiffness!(K, dh, cv, C)
    if bottom_bc != :rough_platens
        assemble_traction!(f, dh, fv, top, Vec(0.0, 0.0, -σ0))
    end

    # --- solve ---------------------------------------------------------------
    f_ext = copy(f)                # apply! rewrites f for inhomogeneous Dirichlet BCs
    apply!(K, f, ch)
    u = K \ f
    apply!(u, ch)

    # residual r = K u - f_ext gives the reaction force at constrained dofs.
    # Separate constraint handlers identify the bottom-face and top-face dofs;
    # dofs are numbered in (u,v,w) triples, so every third dof is a w-component.
    Kfull = allocate_matrix(dh); assemble_stiffness!(Kfull, dh, cv, C)
    ch_bottom = ConstraintHandler(dh); add!(ch_bottom, Dirichlet(:u, bottom, (x, t) -> 0.0, [3])); close!(ch_bottom)
    ch_top    = ConstraintHandler(dh); add!(ch_top,    Dirichlet(:u, top,    (x, t) -> 0.0, [3])); close!(ch_top)
    reaction_z(u) = begin
        r = Kfull * u - f_ext
        Rb = sum(r[d] for d in ch_bottom.prescribed_dofs)
        Rt = sum(r[d] for d in ch_top.prescribed_dofs)
        return Rb, Rt
    end

    if bottom_bc == :rough_platens
        _, Rt = reaction_z(u)                # force platen exerts on the cube (negative)
        u .*= (σ0 * H^2) / abs(Rt)           # scale so the platen carries σ0 H²
    end

    # --- post-process --------------------------------------------------------
    ε_qp, σ_qp = compute_qp_strain_stress(dh, cv, u, C)
    proj = L2Projector(Lagrange{RefHexahedron,1}(), grid)
    σ_proj  = project(proj, σ_qp, qr)                    # ordered by the projector's dofs
    ε_proj  = project(proj, ε_qp, qr)
    σ_nodes = evaluate_at_grid_nodes(proj, σ_proj)       # re-ordered by grid node number
    ε_nodes = evaluate_at_grid_nodes(proj, ε_proj)

    # nodal displacement field as a Vector{Vec{3}}
    u_nodes = evaluate_at_grid_nodes(dh, u, :u)

    # --- summary numbers -----------------------------------------------------
    nodes = getnodes(grid)
    top_ids  = [i for (i, n) in enumerate(nodes) if abs(n.x[3] - H) < 1e-8]

    w_top  = [u_nodes[i][3] for i in top_ids]

    # width change at mid-height, measured at the face centres (frame independent)
    node_at(p) = findfirst(n -> norm(n.x - p) < 1e-8, nodes)
    xL, xR = node_at(Vec(0.0, H/2, H/2)), node_at(Vec(H, H/2, H/2))
    yF, yB = node_at(Vec(H/2, 0.0, H/2)), node_at(Vec(H/2, H, H/2))
    ΔW_x = u_nodes[xR][1] - u_nodes[xL][1]
    ΔW_y = u_nodes[yB][2] - u_nodes[yF][2]

    Δw_exact     = σ0 * H / E            # uniaxial shortening
    ΔW_exact     = ν * σ0 * H / E        # uniaxial increase in width

    Rz, Rz_top = reaction_z(u)              # bottom reaction should equal σ0 H²

    if verbose
        println("=" ^ 72)
        println("Bottom BC: $bottom_bc   mesh: $(nel)^3 hexahedra, order $order, $(ndofs(dh)) dofs")
        println("-" ^ 72)
        @printf("Exact uniaxial shortening  Δw = σ0 H / E          = %.4e mm\n", Δw_exact)
        @printf("FE top-face shortening     mean = %.4e  min = %.4e  max = %.4e mm\n",
                -sum(w_top) / length(w_top), -maximum(w_top), -minimum(w_top))
        @printf("Exact uniaxial width increase  ΔW = ν σ0 H / E  = %.4e mm\n", ΔW_exact)
        @printf("FE width increase at mid-height   x-dir = %.4e   y-dir = %.4e mm\n", ΔW_x, ΔW_y)
        @printf("Applied load  σ0 H² = W_can  = %.4f N\n", σ0 * H^2)
        @printf("Vertical reaction at bottom Σ Rz = %.4f N\n", Rz)
        bottom_bc == :rough_platens && @printf("Vertical force from top platen   = %.4f N\n", Rz_top)
        szz = [s[3, 3] for s in σ_nodes]
        @printf("Nodal σzz range: %.4f … %.4f kPa  (applied -%.4f kPa)\n",
                1e3 * minimum(szz), 1e3 * maximum(szz), 1e3 * σ0)
        println()
    end

    # --- VTK output (open in ParaView, use "Warp By Vector" on u) -------------
    fname = "cube_$(bottom_bc)"
    VTKGridFile(fname, dh) do vtk
        write_solution(vtk, dh, u)
        write_projection(vtk, proj, σ_proj, "sigma")
        write_projection(vtk, proj, ε_proj, "epsilon")
        Ferrite.write_constraints(vtk, ch)
    end
    verbose && println("wrote $fname.vtu\n")

    return (; grid, dh, u, u_nodes, σ_nodes, ε_nodes, Δw_exact, ΔW_exact, w_top, ΔW_x, ΔW_y)
end

# ---------------------------------------------------------------------------
# Deformed-edge profile: lateral displacement of the vertical edge at (x=H, y=H/2)
# ---------------------------------------------------------------------------
function edge_profile(res)
    nodes = getnodes(res.grid)
    ids = [i for (i, n) in enumerate(nodes) if abs(n.x[1] - H) < 1e-8 && abs(n.x[2] - H / 2) < 1e-8]
    sort!(ids, by = i -> nodes[i].x[3])
    return [(nodes[i].x[3], res.u_nodes[i][1], res.u_nodes[i][3]) for i in ids]
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    @printf("\nThermoplastic-rubber cube  H = %.0f mm,  E = %.1f MPa,  ν = %.2f\n", H, E, ν)
    @printf("Load: 5 oz tuna can, m = %.4f kg, W = %.4f N  ->  σ0 = W / H² = %.4e MPa = %.3f kPa\n\n",
            m_can, W_can, σ0, 1e3 * σ0)

    res_free  = solve_cube(:frictionless)
    res_clamp = solve_cube(:clamped)
    res_plat  = solve_cube(:rough_platens)

    println("Deformed profile of the side face centre-line (x = H, y = H/2):")
    println("   z [mm]   u_free [mm]   w_free [mm]   u_clamped [mm]  w_clamped [mm]  u_platens [mm]  w_platens [mm]")
    for ((z, uf, wf), (_, uc, wc), (_, up, wp)) in zip(edge_profile(res_free), edge_profile(res_clamp), edge_profile(res_plat))
        @printf("  %6.2f   %10.3e   %10.3e   %12.3e   %12.3e   %13.3e   %13.3e\n", z, uf, wf, uc, wc, up, wp)
    end
end
