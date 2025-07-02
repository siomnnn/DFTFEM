using Gridap
using Gridap.Arrays
using Gridap.CellData
using Gridap.Geometry
using Gridap.Helpers
using ForwardDiff
using NearestNeighbors
using StaticArrays

import Gridap.Arrays: getindex!, evaluate!
import Gridap.CellData: distance

function Gridap.CellData.evaluate!(cache,f::CellField,x::Point)
  cache1,cache2 = cache

  fixed_cache = (KDTreeSearch(;num_nearest_vertices=2, cache1[1].tol), cache1[2:end]...)

  cell_f_cache, f_cache, cell_f, f₀ = cache2
  @check f === f₀ "Wrong cache"

  cell = Gridap.CellData._point_to_cell!(fixed_cache, x)
  cf = getindex!(cell_f_cache, cell_f, cell)
  fx = evaluate!(f_cache, cf, x)
  return fx
end

function Gridap.CellData._point_to_cell_cache(searchmethod::KDTreeSearch,trian::Triangulation)
  model = get_active_model(trian)
  topo = get_grid_topology(model)
  vertex_coordinates = Gridap.Geometry.get_node_coordinates(model)
  kdtree = KDTree(map(nc -> SVector(Tuple(nc)), vertex_coordinates))
  D = num_cell_dims(trian)
  cell_to_vertices = get_cell_node_ids(model)
  vertex_to_cells = Arrays.inverse_table(cell_to_vertices)
  cell_to_ctype = get_cell_type(trian)
  ctype_to_polytope = get_polytopes(trian)
  cell_map = get_cell_map(trian)
  table_cache = array_cache(vertex_to_cells)
  cache1 = searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache
end

#function _point_to_cell!(cache, x::Point)
#  searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache = cache
#
#  function cell_distance(cell::Integer)
#    ctype = cell_to_ctype[cell]
#    polytope = ctype_to_polytope[ctype]
#    cmap = cell_map[cell]
#    inv_cmap = inverse_map(cmap)
#    return distance(polytope, inv_cmap, x)
#  end
#  
#  # Find the nearest vertices to the point `x` in the triangulation
#  vertices, distances = knn(kdtree, get_array(ForwardDiff.value(x)), searchmethod.num_nearest_vertices, true)
#
#  T = eltype(distances)
#  tol = max(1000*eps(T), T(searchmethod.tol))
#  for vertex in vertices
#
#    # Find all neighbouring cells
#    cells = getindex!(table_cache,vertex_to_cells,vertex)
#    @assert !isempty(cells)
#
#    # Calculate the distance from the point to all the neighbor cells. Without
#    # round-off, and with non-overlapping cells, the distance would be
#    # negative for exactly one cell and positive for all other ones. Due
#    # to round-off, the distance can be slightly negative or slightly
#    # positive for points near cell boundaries, in particular near
#    # vertices. In this case, choose the cell with the smallest
#    # distance, and check that the distance (if positive) is at most at
#    # round-off level.
#
#    cell = zero(eltype(cells))
#    dist = T(Inf)
#    for jcell in cells
#      jdist = cell_distance(jcell)
#      if jdist < dist
#        cell = jcell
#        dist = jdist
#      end
#    end
#
#    (dist < tol) && return cell
#  end
#
#  # Output error message if cell not found
#  @check false "Point $x is not inside any active cell"
#end

#function Gridap.CellData._point_to_cell!(cache, x::Point)
#  searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache = cache
#
#  function cell_distance(cell::Integer)
#    ctype = cell_to_ctype[cell]
#    polytope = ctype_to_polytope[ctype]
#    cmap = cell_map[cell]
#    inv_cmap = inverse_map(cmap)
#    return distance(polytope, inv_cmap, x)
#  end
#
#  # Find the nearest vertices to the point `x` in the triangulation
#  vertices, distances = knn(kdtree, get_array(ForwardDiff.value(x)), searchmethod.num_nearest_vertices, true)
#
#  T = eltype(distances)
#  tol = max(1000*eps(T), T(searchmethod.tol))
#  for vertex in vertices
#
#    # Find all neighbouring cells
#    #cells = getindex!(table_cache,vertex_to_cells,vertex)
#    cells = filter(i -> vertex in vertex_to_cells[i], 1:size(vertex_to_cells, 1))
#    @assert !isempty(cells)
#
#    # Calculate the distance from the point to all the neighbor cells. Without
#    # round-off, and with non-overlapping cells, the distance would be
#    # negative for exactly one cell and positive for all other ones. Due
#    # to round-off, the distance can be slightly negative or slightly
#    # positive for points near cell boundaries, in particular near
#    # vertices. In this case, choose the cell with the smallest
#    # distance, and check that the distance (if positive) is at most at
#    # round-off level.
#
#    cell = zero(eltype(cells))
#    dist = T(Inf)
#    for jcell in cells
#      jdist = cell_distance(jcell)
#      if jdist < dist
#        cell = jcell
#        dist = jdist
#      end
#    end
#    #println("final cell: ", cell, ", final dist: ", dist)
#    (dist < tol) && return cell
#  end
#
#  # Output error message if cell not found
#  @check false "Point $x is not inside any active cell"
#end