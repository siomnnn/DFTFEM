using DFTK
using Plots
using Unitful
using UnitfulAtomic
using PseudoPotentialData
using WriteVTK
using LinearAlgebra

# 1. Define lattice and atomic positions
a = 5.431u"angstrom"          # Silicon lattice constant
lattice = a / 2 * [[0 1 1.];  # Silicon lattice vectors
                   [1 0 1.];  # specified column by column
                   [1 1 0.]];

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

model = Model(lattice, DFTK.Element[], AbstractVector[]; n_electrons = 8, terms = [Kinetic(), ExternalFromReal(r -> norm(r .- austrip(a)/2)^2), Hartree(), LDA()])  # Define model with LDA functional and external potential
kgrid = [1, 1, 1]     # k-point grid (Regular Monkhorst-Pack grid)
Ecut = 7              # kinetic energy cutoff
# Ecut = 190.5u"eV"  # Could also use eV or other energy-compatible units
basis = PlaneWaveBasis(model; Ecut, kgrid)

# 3. Run the SCF procedure to obtain the ground state
scfres = self_consistent_field(basis, tol=1e-5);
#println(typeof(scfres.τ))
#densities = DFTK.LibxcDensities(basis, 1, scfres.ρ, scfres.τ)
save_scfres("julia/dftk_results.vts", scfres; save_ψ=true)
scfres.energies

#@NamedTuple{
#    ham::Hamiltonian,
#    basis::PlaneWaveBasis{Float64, Float64, DFTK.CPU, FFTGrid{Float64, Float64, Array{SVector{3, Int64}, 3}, Array{SVector{3, Float64}, 3}}, Vector{SVector{3, Int64}}},
#    energies::Energies{Float64},
#    converged::Bool,
#    occupation_threshold::Float64,
#    ρ::Array{Float64, 4},
#    τ::Nothing,
#    α::Float64,
#    eigenvalues::Vector{Vector{Float64}},
#    occupation::Vector{Vector{Float64}},
#    εF::Float64,
#    n_bands_converge::Int64,
#    n_iter::Int64,
#    n_matvec::Int64,
#    ψ::Vector{Matrix{ComplexF64}},
#    diagonalization::Vector{
#        @NamedTuple{
#            λ::Vector{Vector{Float64}},
#            X::Vector{Matrix{ComplexF64}},
#            residual_norms::Vector{Vector{Float64}},
#            n_iter::Vector{Int64},
#            converged::Bool,
#            n_matvec::Int64
#        }
#    },
#    stage::Symbol,
#    history_Δρ::Vector{Float64},
#    history_Etot::Vector{Float64},
#    timedout::Bool,
#    mixing::χ0Mixing,
#    runtime_ns::UInt64,
#    algorithm::String
#}

#stack(scfres.eigenvalues)
#stack(scfres.occupation)

#rvecs = collect(r_vectors(basis))[:, 1, 1]  # slice along the x axis
#x = [r[1] for r in rvecs]                   # only keep the x coordinate
#plot(x, scfres.ρ[1, :, 1, 1], label="", xlabel="x", ylabel="ρ", marker=2)
#
#compute_forces_cart(scfres)
#
#bands = compute_bands(scfres, MonkhorstPack(6, 6, 6))
#plot_bandstructure(scfres; kline_density=10)
#plot_dos(bands; temperature=1e-3, smearing=Smearing.FermiDirac())