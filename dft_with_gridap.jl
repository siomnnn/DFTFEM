using DFTK
using Unitful
using UnitfulAtomic
using PseudoPotentialData
using Libxc
using DftFunctionals

import Gmsh: gmsh
using Gridap
using Gridap.FESpaces
using GridapGmsh

using Arpack
using SparseArrays

using Plots

include("gridap_bugfix.jl")

function generate_mesh(lattice::AbstractMatrix{T}, N::Int) where T
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal",0)
    unitless_lattice = austrip.(lattice)

    gmsh.model.add("unit_cell")
    target_mesh_size = (minimum(norm, [unitless_lattice[:, i] for i in 1:3]))/N

    offset = (unitless_lattice[:, 1] + unitless_lattice[:, 2] + unitless_lattice[:, 3])/2

    # first face
    gmsh.model.geo.addPoint(([0, 0, 0] - offset)..., target_mesh_size, 1)
    gmsh.model.geo.addPoint((unitless_lattice[:, 1] - offset)..., target_mesh_size, 2)
    gmsh.model.geo.addPoint((unitless_lattice[:, 1] + unitless_lattice[:, 2] - offset)..., target_mesh_size, 3)
    gmsh.model.geo.addPoint((unitless_lattice[:, 2] - offset)..., target_mesh_size, 4)

    gmsh.model.geo.addLine(1, 2, 5)
    gmsh.model.geo.addLine(2, 3, 6)
    gmsh.model.geo.addLine(3, 4, 7)
    gmsh.model.geo.addLine(4, 1, 8)

    gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 9)
    gmsh.model.geo.addPlaneSurface([9], 10)

    # extrude to volume
    gmsh.model.geo.extrude([2, 10], unitless_lattice[:, 3]...)
    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(3, [1], 1, "unit_cell")

    # specify periodicity and generate mesh
    translation = [1, 0, 0, unitless_lattice[1, 3],
                   0, 1, 0, unitless_lattice[2, 3],
                   0, 0, 1, unitless_lattice[3, 3],
                   0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(2, [32], [10], translation)
    translation = [1, 0, 0, unitless_lattice[1, 2],
                   0, 1, 0, unitless_lattice[2, 2],
                   0, 0, 1, unitless_lattice[3, 2],
                   0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(2, [27], [19], translation)
    translation = [1, 0, 0, unitless_lattice[1, 1],
                   0, 1, 0, unitless_lattice[2, 1],
                   0, 0, 1, unitless_lattice[3, 1],
                   0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(2, [23], [31], translation)
    gmsh.model.mesh.generate(3)
    gmsh.write("julia/unit_cell.msh")

    gmsh.finalize()
end

function get_hartree_pot(model, Ω, dΩ, N_el, vol, ρ)
    reffe = ReferenceFE(lagrangian, Float64, 2)
    V = TestFESpace(model, reffe; conformity=:H1, constraint=:zeromean)
    U = TrialFESpace(V)

    a(u, v) = ∫(∇(u) ⋅ ∇(v))dΩ

    mean_free_ρ(x) = ρ(x) - N_el/vol
    l(v) = 4*pi*(∫(mean_free_ρ*v)dΩ)

    op = AffineFEOperator(a, l, V, U)

    ls = BackslashSolver()
    solver = LinearFESolver(ls)
    V_hartree = solve(solver, op)
    
    writevtk(
        Ω,
        "julia/hartree_pot";
        cellfields=["pot" => V_hartree],
    )

    return V_hartree
end

function get_Vxc(xc, ρ)
    functionals = xc.functionals

    function Vxc(x)
        #return -(3/pi)^(1/3)*ρ(x)^(1/3)
        Vρ = 0
        for functional in functionals
            _, V_new = DftFunctionals.potential_terms(functional, ρ(x)*ones(1, 1))
            Vρ += V_new[1, 1]
        end
        return Vρ
    end

    return Vxc
end

function solve_linear_FEM(model, Ω, dΩ, k, V_eff, n_eigvals::Int)
    reffe = ReferenceFE(lagrangian, Float64, 1)
    V = FESpace(model, reffe; conformity=:H1, vector_type=Vector{ComplexF64})
    U = V
    
    total_pot(x) = 0.5(k ⋅ k) + V_eff(x)

    a(u, v) = ∫(0.5*conj(∇(u)) ⋅ ∇(v))dΩ + ∫(0.5im*(v*(k ⋅ conj(∇(u))) - conj(u)*(k ⋅ ∇(v))))dΩ + ∫(total_pot*conj(u)*v)dΩ
    b(u, v) = ∫(conj(u)*v)dΩ

    A = assemble_matrix(a, U, V)
    B = assemble_matrix(b, U, V)

    λ, ϕ = eigs(A, B; nev=n_eigvals, which=:SR)
    return λ, ϕ
end

function solve_linear_ks(model, Ω, dΩ, k, V_ext, N_el::Int, vol, ρ::FEFunction, n_eigvals::Int, xc::DFTK.Xc)
    V_hartree = get_hartree_pot(model, Ω, dΩ, N_el, vol, ρ)
    Vxc = get_Vxc(xc, ρ)
    potential(x) = V_ext(x) + V_hartree(x) + N_el/vol + Vxc(x)
    
    return solve_linear_FEM(model, Ω, dΩ, k, potential, n_eigvals)
end

function next_density(model, Ω, orbitals, N_el)
    reffe = ReferenceFE(lagrangian, Float64, 1)
    U_real = FESpace(model, reffe; conformity=:H1)
    quad = CellQuadrature(Ω, 2)
    orb_norms = [sqrt(sum(integrate(FEFunction(U_real, abs.(orbitals[:, i]).^2), quad))) for i in 1:N_el]
    
    return sum([2*abs.(orbitals[:, i]./orb_norms[i]).^2 for i in 1:Int(N_el/2)])
end

function solve_ks(model, Ω, k, V_ext, N_el, vol, ρ_start, n_eigvals, functionals)
    solver = DFTK.scf_anderson_solver()
    reffe = ReferenceFE(lagrangian, Float64, 1)
    U_real = FESpace(model, reffe; conformity=:H1)

    num_dofs = num_free_dofs(U_real)

    dΩ = Measure(Ω, 2)

    function fixpoint_map(ρin, info)
        (; orbitals, eigenvalues, n_iter, converged, timedout) = info
        n_iter += 1
        println("n_iter: ", n_iter)

        ρfunc = FEFunction(U_real, ρin, Vector{Float64}([]))
        println("N_el: ", sum(integrate(ρfunc, CellQuadrature(Ω, 2))))

        writevtk(
            Ω,
            "julia/densities";
            cellfields=["density$n_iter" => ρfunc],
        )

        eigenvalues, orbitals = solve_linear_ks(model, Ω, dΩ, k, V_ext, N_el, vol, ρfunc, n_eigvals, functionals)
        
        print(size(orbitals))
        next_ρ = next_density(model, Ω, orbitals, N_el)
        println("energy: ", calculate_energy(eigenvalues, N_el))

        return next_ρ, (; orbitals, eigenvalues, n_iter, converged, timedout)
    end

    info_start = (; orbitals=fill(zeros((num_dofs,)), (Int(N_el/2),)), eigenvalues=[], n_iter=0, converged=false, timedout=false)

    ρout, info = solver(fixpoint_map, get_free_dof_values(ρ_start), info_start; maxiter=10)

    return FEFunction(U_real, ρout, Vector{Float64}([])), info.orbitals, info.eigenvalues
end

function calculate_energy(eigvals, N_el)
    return 2*sum(real.(eigvals[1:Int(N_el/2)]))
end
energies = []
for N in [4, 6, 8, 10, 12, 14, 16, 18, 20]
    #N = 8
    a = 5.431u"angstrom"          # Silicon lattice constant
    lattice = a / 2 * [[0 1 1.];  # Silicon lattice vectors
                       [1 0 1.];  # specified column by column
                       [1 1 0.]];
    #lattice = [[1 0 0]; [0 1 0]; [0 0 1.]];

    unitless_lattice = austrip.(lattice)

    N_el = 8

    pd_lda_family = PseudoFamily("dojo.nc.sr.lda.v0_4_1.standard.upf")
    Si = ElementPsp(:Si, pd_lda_family)
    #Si = ElementPsp(:Si, nothing)

    atoms     = [Si, Si]
    positions = [ones(3)/8, -ones(3)/8]
    real_space_positions_au = [unitless_lattice*position for position in positions]

    V_ext_preliminary(r) = sum([DFTK.local_potential_real(atom, norm(VectorValue(r...) - VectorValue(position...))) for (atom, position) in zip(atoms, real_space_positions_au)])

    V_ext(r) = 0.1*norm(r)^2

    k = VectorValue(0, 0, 0)

    n_eigvals = 26

    generate_mesh(lattice, N)
    model = GmshDiscreteModel("julia/unit_cell.msh")
    Ω = Triangulation(model)
    dΩ = Measure(Ω, 2)

    reffe = ReferenceFE(lagrangian, Float64, 1)
    U_real = FESpace(model, reffe; conformity=:H1)

    vol = austrip(DFTK.compute_unit_cell_volume(lattice))
    ρ = interpolate(x -> N_el/vol, U_real)
    xc = LDA()

    Vxc = get_Vxc(xc, ρ)
    points = range(-austrip(a)/2, austrip(a)/2, 100)
    plot(points, Vxc.(VectorValue.(points, points, points)))
    savefig("julia/xc_potential.png")

    ρ, orbitals, eigvals = solve_ks(model, Ω, k, V_ext, N_el, vol, ρ, n_eigvals, xc)

    println("energy: ", calculate_energy(eigvals, N_el))
    global energies = push!(energies, calculate_energy(eigvals, N_el))

    println("done")
end
open("julia/energies.txt", "w") do f
    println(f, "Eigenvalue sums for N = [4, 6, 8, 10, 12, 14, 16, 18, 20]:")
    println(f, "---------------------------------")
    for energy in energies
        println(f, energy)
    end
end
#plot([4, 6, 8, 10, 12, 14, 16], abs.(energies .- 36.38747469769667), xlabel="N", ylabel="Energy (Hartree)", title="Energy vs N", legend=false, xscale=:log10, yscale=:log10)
#savefig("julia/energy_vs_N.png")