using Gridap

N = 100
model = CartesianDiscreteModel((0, π, 0, π), (N, N))

reffe = ReferenceFE(lagrangian, Float64, 1)
V = FESpace(model, reffe; conformity=:H1, dirichlet_tags="boundary")
U = V

x = Triangulation(model)
dx = Measure(x, 2)

a(u, v) = ∫(∇(u) ⋅ ∇(v))dx
b(v) = 0.0

A = AffineFEOperator(a, b, U, V).op.matrix

using Arpack

nev = 10
λ, ϕ = eigs(A; nev=nev, which=:SM)

writevtk(
    x,
    "julia/modes";
    cellfields=["mode$i" => FEFunction(U, real(ϕ[:, i])) for i in 1:nev],
)

println(λ[1:nev])