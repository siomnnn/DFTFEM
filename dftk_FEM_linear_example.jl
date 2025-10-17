# taken straight from the DFTK docs, except different pseudopotential and FEM calculations

using DFTK
using Plots
using Unitful
using UnitfulAtomic
using PseudoPotentialData
using StaticArrays

# 1. Define lattice and atomic positions
a = 5.431u"angstrom"          # Silicon lattice constant
lattice = a / 2 * [[0 1 1.];  # Silicon lattice vectors
                   [1 0 1.];  # specified column by column
                   [1 1 0.]];

family_gth = PseudoFamily("cp2k.nc.sr.lda.v0_1.semicore.gth")
Si = ElementPsp(:Si, family_gth)

# Specify type and positions of atoms
atoms     = [Si, Si]
positions = [ones(3)/8, -ones(3)/8]

# 2. Select model and basis
model = model_DFT(lattice, atoms, positions; functionals=LDA())
kgrid = [1, 1, 1]     # k-point grid (Regular Monkhorst-Pack grid)
Ecut = 100            # kinetic energy cutoff
# Ecut = 190.5u"eV"  # Could also use eV or other energy-compatible units
basis = PlaneWaveBasis(model; Ecut, kgrid)#, kshift=[1/4, 1/4, 1/4], use_symmetries_for_kpoint_reduction=false)
# Note the implicit passing of keyword arguments here:
# this is equivalent to PlaneWaveBasis(model; Ecut=Ecut, kgrid=kgrid)


ρ_function = v -> (sin((2π*basis.model.inv_lattice*v)[1])^2 * sin((2π*basis.model.inv_lattice*v)[2])^2 * sin((2π*basis.model.inv_lattice*v)[3])^2 * 32/basis.model.unit_cell_volume + 4/basis.model.unit_cell_volume)

ψ_pw = nothing#DFTK.random_orbitals(basis, basis.kpoints[1], 7)
println(size(r_vectors(basis)))
println(size(r_vectors_cart(basis)))
ρ_pw = reshape(ρ_function.(r_vectors_cart(basis)), (size(r_vectors_cart(basis))..., 1))
println(size(ρ_pw))
E_pw, ham_pw = DFTK.energy_hamiltonian(basis, ψ_pw, [[1.0, 1.0, 1.0, 1.0, 0.0, 0.0]]; ρ=ρ_pw)
diag_res = DFTK.next_density(ham_pw)
println("diag eigenvalues: ", diag_res.eigenvalues)

println("starting FEM calculation")
h = a/40

DFTK.reset_timer!(DFTK.timer)
basis_fem = FiniteElementBasis(model; h, kgrid=[1, 1, 1])#, kshift=[1/4, 1/4, 1/4])
println("finished basis setup")

ψ_fem = nothing#DFTK.random_orbitals(basis_fem, basis.kpoints[1], 7)
ρ_fem = reshape(ρ_function.(DFTK.get_free_dof_positions(basis_fem, :ρ)), (DFTK.get_n_free_dofs(basis_fem, :ρ), 1))
constraint_matrix = DFTK.get_constraint_matrix(basis_fem, :ρ)
density_integral = (real(ρ_fem' * constraint_matrix') * DFTK.get_overlap_matrix(basis_fem, :ρ) * real(constraint_matrix * ones(size(ρ_fem, 1))))[1]
println("integral of density: ", density_integral)
E_fem, ham_fem = DFTK.energy_hamiltonian(basis_fem, ψ_fem, [[1.0, 1.0, 1.0, 1.0, 0.0, 0.0]]; ρ=ρ_fem/(density_integral/8))
diag_res = DFTK.next_density(ham_fem)
println("diag eigenvalues: ", diag_res.eigenvalues)

@show DFTK.timer