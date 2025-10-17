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
kgrid = [2, 2, 2]     # k-point grid (Regular Monkhorst-Pack grid)
Ecut = 100            # kinetic energy cutoff
# Ecut = 190.5u"eV"  # Could also use eV or other energy-compatible units
basis = PlaneWaveBasis(model; Ecut, kgrid)#, kshift=[1/4, 1/4, 1/4], use_symmetries_for_kpoint_reduction=false)
# Note the implicit passing of keyword arguments here:
# this is equivalent to PlaneWaveBasis(model; Ecut=Ecut, kgrid=kgrid)

# 3. Run the SCF procedure to obtain the ground state
scfres = self_consistent_field(basis, tol=1e-5)
println("eigenvalues: ", scfres.eigenvalues)
println("energy breakdown from planewave reference: ", DFTK.todict(scfres.energies))

println("starting FEM calculation")
h = a/40
degree = 1

DFTK.reset_timer!(DFTK.timer)
basis_fem = basis = FiniteElementBasis(model; h, kgrid=[1, 1, 1])#, kshift=[1/4, 1/4, 1/4])
println("finished basis setup")

scfres_fem = self_consistent_field(basis_fem, tol=1e-5)
println("eigenvalues: ", scfres_fem.eigenvalues)
println("energy breakdown from FEM: ", DFTK.todict(scfres_fem.energies))
@show DFTK.timer