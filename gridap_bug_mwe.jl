using Gridap
using Gridap.CellData
using Gridap.Geometry
using GridapGmsh

import Gridap.Geometry: get_node_coordinates

include("gridap_bugfix.jl")

#function Gridap.CellData._point_to_cell_cache(searchmethod::KDTreeSearch,trian::Triangulation)
#  model = get_active_model(trian)
#  topo = get_grid_topology(model)
#  vertex_coordinates = get_node_coordinates(model)
#  kdtree = KDTree(map(nc -> SVector(Tuple(nc)), vertex_coordinates))
#  D = num_cell_dims(trian)
#  vertex_to_cells = get_cell_node_ids(model)
#  cell_to_ctype = get_cell_type(trian)
#  ctype_to_polytope = get_polytopes(trian)
#  cell_map = get_cell_map(trian)
#  table_cache = array_cache(vertex_to_cells)
#  cache1 = searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache
#end

model = GmshDiscreteModel("periodic_mesh.msh")
Ω = Triangulation(model)
#dΩ = Measure(Ω, 2)
#
#reffe = ReferenceFE(lagrangian, Float64, 1)
#V = TestFESpace(model, reffe; conformity=:H1, constraint=:zeromean)
#U = TrialFESpace(V)
#
#a(u, v) = ∫(∇(u) ⋅ ∇(v))dΩ
#l(v) = 0
#
#op = AffineFEOperator(a, l, V, U)
#
#ls = LUSolver()
#solver = LinearFESolver(ls)
#solution = solve(solver, op)
#
#solution(VectorValue(0.38983887647284665, -0.15055445742314805, 0.2548498884962602))

field = CellField(x -> 1, Ω)
print("YES WE FINISHED!!!!!   ", evaluate_fixed(field, VectorValue(0.38983887647284665, -0.15055445742314805, 0.2548498884962602)))
#print("YES WE FINISHED!!!!!", field(VectorValue(0.38983887647284665, -0.15055445742314805, 0.2548498884962602)))