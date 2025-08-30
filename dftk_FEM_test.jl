using DFTK
using Plots
using Unitful
using UnitfulAtomic
using PseudoPotentialData
using WriteVTK
using LinearAlgebra
using Ferrite
using FerriteGmsh
using Arpack
using StaticArrays
using Krylov
using SparseArrays
using JLD2

#RECOMPILE DFTK USING Pkg.develop(PackageSpec(path="/home/siomn/nc-simon/uni/bachelorarbeit/DFTK.jl"))

println("starting...")

# 1. Define lattice and atomic positions
a = 5.431u"angstrom"          # Silicon lattice constant
#lattice = a / 2 * [[0 1 1.];  # Silicon lattice vectors
#                   [1 0 1.];  # specified column by column
#                   [1 1 0.]];
lattice = a * [[1 0 0.];
               [0 1 0.];
               [0 0 1.]];

pd_lda_family = PseudoFamily("dojo.nc.sr.lda.v0_4_1.standard.upf")
Si = ElementPsp(:Si, pd_lda_family)

# Specify type and positions of atoms
atoms     = [Si, Si]
positions = [ones(3)/8, -ones(3)/8]

## 2. Select model and basis
#model = model_DFT(lattice, atoms, positions; functionals=LDA())
#kgrid = [1, 1, 1]     # k-point grid (Regular Monkhorst-Pack grid)
#Ecut = 7              # kinetic energy cutoff
## Ecut = 190.5u"eV"  # Could also use eV or other energy-compatible units
#basis = PlaneWaveBasis(model; Ecut, kgrid)
## Note the implicit passing of keyword arguments here:
## this is equivalent to PlaneWaveBasis(model; Ecut=Ecut, kgrid=kgrid)

model = Model(lattice, DFTK.Element[], AbstractVector[]; n_electrons = 8, terms = [Kinetic(), ExternalFromReal(r -> 0.01norm(r .- austrip(a)/2)^2), Hartree(), LDA()])  # Define model with LDA functional and external potential
h = a/20
#h = 1/20

#grid = load_grid_from_file("mesh.msh")
basis = FiniteElementBasis(model; h, degree=1, kpoints=[DFTK.FEMKpoint(1, SVector(0, 0, 0.))], kweights=[1.], precompute_laplacian=true)
println("finished basis setup")

function test_scf()
    pw_basis = PlaneWaveBasis(model; Ecut=7, kgrid=[1, 1, 1])

    println("starting planewave reference SCF...")
    pw_scfres = self_consistent_field(pw_basis)

    println("energy breakdown from planewave reference: ", DFTK.todict(pw_scfres.energies))

    println("starting FEM SCF...")
    scfres = self_consistent_field(basis)

    ρ = scfres.ρ

    println("energy breakdown from FEM: ", DFTK.todict(scfres.energies))

    #out = zeros(eltype(basis), DFTK.get_n_dofs(basis, :ρ))
    #out[DFTK.get_free_dofs(basis, :ρ)] .= DFTK.total_density_FEM(ρ)
    #Ferrite.apply!(out, DFTK.get_constraint_handler(basis, :ρ))
    #VTKGridFile("density", DFTK.get_dof_handler(basis, :ρ)) do vtk
    #    write_solution(vtk, DFTK.get_dof_handler(basis, :ρ), out)
    #end
    open("scfres.txt", "w") do f
        println(f, "iteration, energy history, Δρ history:")
        println(f, "---------------------------------")
        for (i, energy) in enumerate(scfres.history_Etot)
            println(f, "n = $i: ", energy, ", ", scfres.history_Δρ[i])
        end
        println(f, "total runtime: ", scfres.runtime_ns)
        println(f, "number of iterations: ", scfres.n_iter)
        println(f, "energy breakdown: ", DFTK.todict(scfres.energies))
    end
end

test_scf()

function test_hamiltonians()
    println([typeof(term) for term in basis.terms])

    laplace_op = DFTK.NegHalfLaplaceFEMOperator(basis, basis.kpoints[1])
    laplace_op_k = DFTK.NegHalfLaplaceFEMOperator(basis, basis.kpoints[2])

    laplace_eig = DFTK.LOBPCG(laplace_op, zeros(Complex{eltype(basis)}, DFTK.get_n_free_dofs(basis, :ψ), 15), DFTK.get_overlap_matrix(basis, :ψ))
    laplace_eig_k = DFTK.LOBPCG(laplace_op_k, zeros(Complex{eltype(basis)}, DFTK.get_n_free_dofs(basis, :ψ), 15), DFTK.get_overlap_matrix(basis, :ψ))
    eigenvectors = laplace_eig.X
    eigenvectors_k = laplace_eig_k.X
    ψ = [eigenvectors[:, [2, 8, 4]], eigenvectors_k[:, [2, 8, 4]], eigenvectors_k[:, [2, 8, 4]]]
    occupation = [[1, 1, 0], [1, 1, 0], [1, 1, 0]]

    println("starting energy computation...")
    result = energy_hamiltonian(basis, ψ, occupation; ρ=compute_density(basis, ψ, occupation))
    println("energies: ", DFTK.todict(result.energies))
    println("Hamiltonian: ", typeof(result.ham))
    println("direct energies: ", DFTK.todict(DFTK.energy(basis, ψ, occupation; ρ=compute_density(basis, ψ, occupation)).energies))

    println("operators: ", [typeof(op) for op in result.ham.blocks[1].optimized_operators])

    Hψ = result.ham * ψ
    println("size of Hψ: ", size(Hψ))
    println("quadratic form of Hamiltonian: ", sum([dot(ψ[ik][:, n], Hψ[ik][:, n]) * occupation[ik][n] * basis.kweights[ik] for n in eachindex(occupation[1]) for ik in eachindex(ψ)]))

    println("size of total potential: ", size(DFTK.total_local_potential(result.ham)))

    new_H = DFTK.hamiltonian_with_total_potential(result.ham, DFTK.total_local_potential(result.ham))
    println("Hamiltonian acts the same? ", Hψ ≈ new_H * ψ)
end

function test_terms()
    kin = Kinetic()
    term_kin = kin(basis)
    println("kinetic term: ", term_kin)

    laplace_op = DFTK.NegHalfLaplaceFEMOperator(basis, basis.kpoints[1])
    laplace_op_k = DFTK.NegHalfLaplaceFEMOperator(basis, basis.kpoints[2])

    free_dofs = DFTK.get_constraint_handler(basis, :ψ).free_dofs
    println("free dofs: ", size(free_dofs))
    println("size of inverse constraint map: ", size(DFTK.get_inverse_constraint_map(basis, :ψ)))
    println("max free degree in inverse constraint map: ", maximum(DFTK.get_inverse_constraint_map(basis, :ψ)))

    constr = DFTK.get_constraint_matrix(basis, basis.kpoints[1], :ψ)
    constr_k = DFTK.get_constraint_matrix(basis, basis.kpoints[2], :ψ)

    laplace_eig = DFTK.LOBPCG(laplace_op, zeros(Complex{eltype(basis)}, DFTK.get_n_free_dofs(basis, :ψ), 15), DFTK.get_overlap_matrix(basis, :ψ))
    laplace_eig_k = DFTK.LOBPCG(laplace_op_k, zeros(Complex{eltype(basis)}, DFTK.get_n_free_dofs(basis, :ψ), 15), DFTK.get_overlap_matrix(basis, :ψ))
    eigenvectors = laplace_eig.X
    eigenvectors_k = laplace_eig_k.X
    println("eigenvalues of laplace operator at Γ: ", laplace_eig.λ)
    println("lowest eigenvalue should be: ", 0.0)
    println("eigenvalues of laplace operator at L: ", laplace_eig_k.λ)
    println("lowest eigenvalue should be: ", 3/8)
    ψ = [eigenvectors[:, [2, 8, 4]], eigenvectors_k[:, [2, 8, 4]]]
    println("size of ψ: ", size(ψ), ", n_dofs: ", DFTK.get_n_dofs(basis, :ψ))
    occupation = [[1, 1, 0], [1, 1, 0]]
    Eops = DFTK.ene_ops(term_kin, basis, ψ, occupation)
    println(Eops.E)
    println(typeof(Eops.ops[1]))

    ext_direct = DFTK.ExternalFromValues(ones(Complex{eltype(basis)}, DFTK.get_n_free_dofs(basis, :ρ)))
    ext_term_direct = ext_direct(basis)
    println("external term from array: ", typeof(ext_term_direct))

    ext_real = ExternalFromReal(r -> 0.1*norm(r .- austrip(a)/2)^2)
    ext_term_real = ext_real(basis)
    println("external term from function: ", typeof(ext_term_real))

    ρ = compute_density(basis, ψ, occupation)
    println("size of ρ: ", size(ρ))
    println("density correct at ψ dofs? ", sum(abs.(DFTK.reduce_dofs(basis, ρ) - (abs.(ψ[1][:, 1]).^2 + abs.(ψ[1][:, 2]).^2)) .< 1e-8))

    Eops_ext_direct = DFTK.ene_ops(ext_term_direct, basis, ψ, occupation; ρ=ρ)
    println("energy from external term (direct): ", Eops_ext_direct.E)
    println("energy from external term (direct) from operator: ", sum([dot(ψ[1][:, i], Eops_ext_direct.ops[1]*ψ[1][:, i]) * occupation[1][i] for i in 1:size(ψ[1], 2)]))
    println("energy from external term (direct) from operator (nonzero k): ", sum([dot(ψ[2][:, i], Eops_ext_direct.ops[2]*ψ[2][:, i]) * occupation[2][i] for i in 1:size(ψ[2], 2)]))
    ext_constant = ExternalFromReal(r -> 1.0)
    ext_term_constant = ext_constant(basis)
    Eops_ext_constant = DFTK.ene_ops(ext_term_constant, basis, ψ, occupation; ρ=ρ)
    println("energy from external term (constant 1): ", Eops_ext_constant.E)
    println("energy from external term (constant 1) from operator: ", sum([dot(ψ[1][:, i], Eops_ext_constant.ops[1]*ψ[1][:, i]) * occupation[1][i] for i in 1:size(ψ[1], 2)]))
    println("energy from external term (constant 1) from operator (nonzero k): ", sum([dot(ψ[2][:, i], Eops_ext_constant.ops[2]*ψ[2][:, i]) * occupation[2][i] for i in 1:size(ψ[2], 2)]))
    Eops_ext_real = DFTK.ene_ops(ext_term_real, basis, ψ, occupation; ρ=ρ)
    println("energy from external term (real): ", Eops_ext_real.E)
    println("energy from external term (real) from operator: ", sum([dot(ψ[1][:, i], Eops_ext_real.ops[1]*ψ[1][:, i]) * occupation[1][i] for i in 1:size(ψ[1], 2)]))
    println("energy from external term (real) from operator (nonzero k): ", sum([dot(ψ[2][:, i], Eops_ext_real.ops[2]*ψ[2][:, i]) * occupation[2][i] for i in 1:size(ψ[2], 2)]))
    #println("energy from external term (real) from matrix: ", dot(ψ[1][:, 1], (DFTK.Matrix(Eops_ext_real.ops[1]))*ψ[1][:, 1]) * occupation[1][1] + dot(ψ[1][:, 2], (DFTK.Matrix(Eops_ext_real.ops[1]))*ψ[1][:, 2]) * occupation[1][2])
    Eops_ones = DFTK.ene_ops(ext_term_real, basis, ψ, occupation; ρ=reshape(ones(eltype(basis), DFTK.get_n_free_dofs(basis, :ρ)), DFTK.get_n_free_dofs(basis, :ρ), 1))
    println("energy from external term (real with ρ=ones): ", Eops_ones.E)
    println("energy from external term (real with ρ=ones) from operator: ", dot(ones(eltype(basis), DFTK.get_n_free_dofs(basis, :ψ)), (Eops_ones.ops[1]*ones(eltype(basis), DFTK.get_n_free_dofs(basis, :ψ)))))
    println("energy from external term (real with ρ=ones) from operator (nonzero k): ", dot(ones(eltype(basis), DFTK.get_n_free_dofs(basis, :ψ)), (Eops_ones.ops[2]*ones(eltype(basis), DFTK.get_n_free_dofs(basis, :ψ)))))

    #potential = ext_term_real.potential_values
    #println("size of potential: ", size(potential))
    #println("maximum of potential: ", maximum(abs.(potential)))
    #println("number of maxs: ", sum(abs.(potential) .> 3/4 - 1e-15))
    #dof_coords = DFTK.get_dof_positions(basis, :ρ)
    #println("size of dof coordinates: ", size(dof_coords))
    #dh = DFTK.get_dof_handler(basis, :ψ)
    #ip = Ferrite.getfieldinterpolation(dh, Ferrite.find_field(dh, :ψ))
    #ref_coords = hcat(Ferrite.reference_coordinates(ip)...)
    #cell_values = DFTK.get_cell_values(basis, :ψ)
    #f = r -> norm(r .- 1/2)^2
    #f_evals = f.(dof_coords)
    #zero_coords = findall(abs.(f_evals) .> 3/4 - 1e-15)
    #println("dofs with max potential: ", zero_coords)
    #println("coords of dofs with max potential: ", dof_coords[zero_coords])
    #println("corresponding real dofs: ", DFTK.get_free_dofs(basis, :ρ)[zero_coords])

    #ρ = ones(eltype(basis), DFTK.get_n_free_dofs(basis, :ρ))
    constraint_matrix = DFTK.get_constraint_matrix(basis, basis.kpoints[1], :ρ)
    (x, stats) = minres_qlp(constraint_matrix' * DFTK.get_neg_half_laplace_matrix(basis, :ρ) * constraint_matrix, DFTK.total_density_FEM(ρ) .- sum(ρ)/DFTK.get_n_free_dofs(basis, :ρ))
    #println("minres_qlp result: ", x)
    println("minres_qlp result size: ", size(x))
    println("minres_qlp result norm: ", norm(x))
    println("minres_qlp result mean: ", sum(x)/DFTK.get_n_free_dofs(basis, :ρ))
    println("minres_qlp stats: ", stats)
    println("minres_qlp converged? ", stats.solved)

    println("mean value comparison: ", basis.model.n_electrons / basis.model.unit_cell_volume, " vs. ", sum(ρ)/DFTK.get_n_free_dofs(basis, :ρ))
    hartree = Hartree()
    term_hartree = hartree(basis)
    println("hartree term: ", typeof(term_hartree))
    Eops_hartree = DFTK.ene_ops(term_hartree, basis, ψ, occupation; ρ=ρ)
    println("energy from hartree term: ", Eops_hartree.E)

    out = zeros(eltype(basis), DFTK.get_n_dofs(basis, :ρ))
    out[DFTK.get_free_dofs(basis, :ρ)] .= DFTK.total_density_FEM(ρ)
    Ferrite.apply!(out, DFTK.get_constraint_handler(basis, :ρ))
    VTKGridFile("density", DFTK.get_dof_handler(basis, :ρ)) do vtk
        write_solution(vtk, DFTK.get_dof_handler(basis, :ρ), out)
    end
    ρ_compare = abs.(ψ[1][:, 1]).^2 + abs.(ψ[1][:, 2]).^2
    out = zeros(eltype(basis), DFTK.get_n_dofs(basis, :ψ))
    out[DFTK.get_free_dofs(basis, :ψ)] .= ρ_compare
    Ferrite.apply!(out, DFTK.get_constraint_handler(basis, :ψ))

    pot = zeros(eltype(basis), DFTK.get_n_dofs(basis, :ρ))
    pot[DFTK.get_free_dofs(basis, :ρ)] .= x
    Ferrite.apply!(pot, DFTK.get_constraint_handler(basis, :ρ))
    VTKGridFile("hartree_pot", DFTK.get_dof_handler(basis, :ρ)) do vtk
        write_solution(vtk, DFTK.get_dof_handler(basis, :ρ), pot)
    end

    xc = LDA()
    term_xc = xc(basis)
    println("xc term: ", typeof(term_xc))
    Eops_xc = DFTK.ene_ops(term_xc, basis, ψ, occupation; ρ=ρ)
    println("energy from xc term: ", Eops_xc.E)

    constr = DFTK.get_constraint_matrix(basis, basis.kpoints[1], :ψ)
    println("check that ψ is not real-valued: ", maximum(abs.(imag.(ψ[2][:, 1]))))
    println("check that R is not real-valued: ", maximum(abs.(imag.(constr))))
    VTKGridFile("density_from_psi", DFTK.get_dof_handler(basis, :ψ)) do vtk
        write_solution(vtk, DFTK.get_dof_handler(basis, :ψ), real.((abs.(constr * ψ[2][:, 1]).^2 + abs.(constr * ψ[2][:, 2]).^2)))
    end
    #VTKGridFile("eigenfunction", DFTK.get_dof_handler(basis, :ψ)) do vtk
    #    write_solution(vtk, DFTK.get_dof_handler(basis, :ψ), abs.(ψ[:, 3]).^2)
    #end
end

function test_ops()
    out = zeros(eltype(basis), n_dofs)

    #noop = DFTK.NoopFEMOperator(basis)
    #DFTK.apply!(out, noop, ones(eltype(basis), n_dofs))
    #H = DFTK.Matrix(noop)
    #println("noop size: ", size(H), ", n_dofs: ", n_dofs)

    multiplication = DFTK.FEMRealSpaceMultiplication(basis, ones(eltype(basis), n_dofs))
    ψ = ones(eltype(basis), n_dofs)
    DFTK.remove_bc!(ψ, DFTK.getconstrainthandler(basis))
    DFTK.apply!(out, multiplication, ψ)
    H = DFTK.Matrix(multiplication)
    println("mult size: ", size(H), ", n_dofs: ", n_dofs)
    out_test = H * ψ
    #Ferrite.apply!(out_test, DFTK.getconstrainthandler(basis))

    out_bc = zeros(eltype(basis), n_dofs)
    ψ = ones(eltype(basis), n_dofs)
    Ferrite.apply!(ψ, DFTK.getconstrainthandler(basis))
    DFTK.apply!(out_bc, multiplication, ψ)

    println("overlap consistent with multiplication by 1? ", maximum(abs.(basis.overlap_matrix - H)) < 1e-14)
    println("application consistent with matrix multiplication after explicit bc? ", all(abs.(out - out_test) .< 1e-12))
    println("max diff: ", maximum(abs.(out_test - out)))
    mult2 = DFTK.FEMRealSpaceMultiplication(basis, 1:n_dofs)
    #println("eigenvalues of multiplication with 1:n_dofs: ", DFTK.LOBPCG(mult2, zeros(eltype(basis), n_dofs, 15), basis.overlap_matrix).λ)

    out = zeros(eltype(basis), n_dofs)

    kinetic = DFTK.NegHalfLaplaceFEMOperator(basis)

    ψ = ones(eltype(basis), n_dofs)
    DFTK.remove_bc!(ψ, DFTK.getconstrainthandler(basis))
    #Ferrite.apply!(ψ, DFTK.getconstrainthandler(basis))
    DFTK.apply!(out, kinetic, ψ)
    #H = DFTK.Matrix(kinetic)
    #println("laplace size: ", size(H), ", n_dofs: ", n_dofs)
    println("result correctly zero? ", all(abs.(out) .< 1e-12))
    println("kinetic energy: ", ψ' * out)

    laplace_eig = DFTK.LOBPCG(basis.neg_half_laplacian, zeros(eltype(basis), n_dofs, 15), basis.overlap_matrix)
    ψ = laplace_eig.X[:, 2]
    #Ferrite.apply!(ψ, DFTK.getconstrainthandler(basis))
    #println("norm of eigenvector with bc: ", norm(ψ, basis))
    #DFTK.apply_bc!(ψ, DFTK.getconstrainthandler(basis))
    #println("norm of eigenvector without bc: ", norm(ψ, basis))
    println("same before and after boundary conditions? ", sum(abs.(ψ - laplace_eig.X[:, 2]) .< 1e-12))
    out_mat = zeros(eltype(basis), n_dofs)
    DFTK.apply!(out_mat, kinetic, ψ)
    DFTK.remove_bc!(out_mat, DFTK.getconstrainthandler(basis))
    #out = basis.neg_half_laplacian * ψ
    println("correct eigenvalues: ", laplace_eig.λ)
    println("kinetic energy: ", ψ' * out_mat)

    no_pregen_basis = FiniteElementBasis(model; h, degree=1, precompute_laplacian=false)
    matrix_free_kinetic = DFTK.NegHalfLaplaceFEMOperator(no_pregen_basis)
    #Ferrite.apply!(ψ, DFTK.getconstrainthandler(no_pregen_basis))
    out_matrix_free = zeros(eltype(no_pregen_basis), DFTK.getndofs(no_pregen_basis))
    DFTK.apply!(out_matrix_free, matrix_free_kinetic, ψ)
    DFTK.remove_bc!(out_matrix_free, DFTK.getconstrainthandler(basis))
    true_result = laplace_eig.λ[8]*(basis.overlap_matrix*ψ)
    DFTK.remove_bc!(true_result, DFTK.getconstrainthandler(basis))
    matmul = basis.neg_half_laplacian*laplace_eig.X[:, 8]
    DFTK.remove_bc!(matmul, DFTK.getconstrainthandler(basis))
    println("matrix free consistent with matrix mult? ", all((out_matrix_free - out_mat) .< 1e-12))
    println("matrix free kinetic energy: ", ψ' * out_matrix_free)
    println("LOBPCG error: ", maximum(abs.(matmul - true_result)))
    println("error for matrix-based: ", maximum(abs.(out_mat - true_result)))
    println("error for matrix-free: ", maximum(abs.(out_matrix_free - true_result)))
    #println("correct eigenvector? ", (out - laplace_eig.λ[8]*(basis.overlap_matrix*ψ)))

    matrix_free_eig = DFTK.LOBPCG(matrix_free_kinetic, zeros(eltype(no_pregen_basis), DFTK.getndofs(no_pregen_basis), 15), basis.overlap_matrix)
    println("matrix free eigenvalues: ", matrix_free_eig.λ)
    #println("matrix free eigenvector:", matrix_free_eig.X[:, 1])
    #println("matrix free kinetic energy using matrix mult: ", matrix_free_eig.X[:, 3]' * (basis.neg_half_laplacian * matrix_free_eig.X[:, 3]))
    #out_ev_test = zeros(eltype(no_pregen_basis), DFTK.getndofs(no_pregen_basis))
    #DFTK.apply!(out_ev_test, matrix_free_kinetic, matrix_free_eig.X[:, 3])
    #println("matrix free kinetic energy using matrix free: ", matrix_free_eig.X[:, 3]' * out_ev_test)

    VTKGridFile("eigenfunction", DFTK.getdofhandler(basis)) do vtk
        write_solution(vtk, DFTK.getdofhandler(basis), abs.(ψ).^2)
    end
end