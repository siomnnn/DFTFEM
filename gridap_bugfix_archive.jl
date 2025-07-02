using Gridap
using Gridap.Arrays
using Gridap.CellData
using Gridap.Helpers
using ForwardDiff
using NearestNeighbors
using StaticArrays

import Gridap.Arrays: getindex!, evaluate!
import Gridap.CellData: _point_to_cell!, distance

function evaluate_fixed(f,x...)
  c = return_cache_fixed(f,x...)
  y = evaluate_fixed!(c,f,x...)
end

return_cache_fixed(f::CellField,x::Point) = return_cache_fixed(Interpolable(f),x)

function return_cache_fixed(a::Interpolable,x::Point)
  f = a.uh
  trian = get_triangulation(f)
  cache1 = _point_to_cell_cache_fixed(a.searchmethod,trian)

  cell_f = get_array(f)
  cell_f_cache = array_cache(cell_f)
  cf = testitem(cell_f)
  f_cache = return_cache(cf,x)
  cache2 = cell_f_cache, f_cache, cell_f, f

  return cache1, cache2
end

function _point_to_cell_cache_fixed(searchmethod::KDTreeSearch,trian::Triangulation)
  model = get_active_model(trian)
  topo = Gridap.Geometry.get_grid_topology(model)
  vertex_coordinates = Gridap.Geometry.get_node_coordinates(model)
  kdtree = KDTree(map(nc -> SVector(Tuple(nc)), vertex_coordinates))
  D = num_cell_dims(trian)
  vertex_to_cells = Gridap.Geometry.get_cell_node_ids(model)
  cell_to_ctype = Gridap.Geometry.get_cell_type(trian)
  ctype_to_polytope = Gridap.Geometry.get_polytopes(trian)
  cell_map = Gridap.Geometry.get_cell_map(trian)
  table_cache = array_cache(vertex_to_cells)
  cache1 = searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache
end

#function _point_to_cell_cache(searchmethod::KDTreeSearch,trian::Triangulation)
#  model = get_active_model(trian)
#  topo = get_grid_topology(model)
#  vertex_coordinates = Geometry.get_vertex_coordinates(topo)
#  kdtree = KDTree(map(nc -> SVector(Tuple(nc)), vertex_coordinates))
#  D = num_cell_dims(trian)
#  vertex_to_cells = Gridap.Geometry.get_faces(topo, 0, D)
#  cell_to_ctype = get_cell_type(trian)
#  ctype_to_polytope = get_polytopes(trian)
#  cell_map = get_cell_map(trian)
#  table_cache = array_cache(vertex_to_cells)
#  cache1 = searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache
#end

function evaluate_fixed!(cache,f::CellField,x::Point)
  cache1,cache2 = cache
  #index = findfirst(item -> item == 216, cache1[2].indices)
  #println(cache1[2].data[index])

  fixed_cache = (KDTreeSearch(;num_nearest_vertices=2, cache1[1].tol), cache1[2:end]...)

  cell_f_cache, f_cache, cell_f, f₀ = cache2
  @check f === f₀ "Wrong cache"

  cell = _point_to_cell_fixed!(fixed_cache, x)
  cf = getindex!(cell_f_cache, cell_f, cell)
  fx = evaluate!(f_cache, cf, x)
  return fx
end

function _point_to_cell_fixed!(cache, x::Point)
  searchmethod, kdtree, vertex_to_cells, cell_to_ctype, ctype_to_polytope, cell_map, table_cache = cache

  function cell_distance(cell::Integer)
    ctype = cell_to_ctype[cell]
    polytope = ctype_to_polytope[ctype]
    cmap = cell_map[cell]
    inv_cmap = inverse_map(cmap)
    return distance(polytope, inv_cmap, x)
  end

  # Find the nearest vertices to the point `x` in the triangulation
  vertices, distances = knn(kdtree, get_array(ForwardDiff.value(x)), searchmethod.num_nearest_vertices, true)
  #println("x: ", x)
  #println("id: ", vertices, ", dist: ", distances)
  #println(searchmethod.num_nearest_vertices)

  T = eltype(distances)
  tol = max(1000*eps(T), T(searchmethod.tol))
  for vertex in vertices

    # Find all neighbouring cells
    #cells = getindex!(table_cache,vertex_to_cells,vertex)
    cells = filter(i -> vertex in vertex_to_cells[i], 1:size(vertex_to_cells, 1))
    #println(cells)
    @assert !isempty(cells)

    # Calculate the distance from the point to all the neighbor cells. Without
    # round-off, and with non-overlapping cells, the distance would be
    # negative for exactly one cell and positive for all other ones. Due
    # to round-off, the distance can be slightly negative or slightly
    # positive for points near cell boundaries, in particular near
    # vertices. In this case, choose the cell with the smallest
    # distance, and check that the distance (if positive) is at most at
    # round-off level.

    cell = zero(eltype(cells))
    dist = T(Inf)
    for jcell in cells
      jdist = cell_distance(jcell)
      if jdist < dist
        cell = jcell
        dist = jdist
      end
    end
    #println("final cell: ", cell, ", final dist: ", dist)
    (dist < tol) && return cell
  end

  # Output error message if cell not found
  @check false "Point $x is not inside any active cell"
end

#Tuple{
#    KDTreeSearch, 
#    NearestNeighbors.KDTree{
#        StaticArraysCore.SVector{3, Float64}, Distances.Euclidean, Float64, StaticArraysCore.SVector{3, Float64}
#    }, 
#    Table{
#        Int32, Vector{Int32}, Vector{Int32}
#    }, 
#    Vector{Int8}, 
#    Vector{Gridap.ReferenceFEs.ExtrusionPolytope{3}}, 
#    LazyArray{
#        FillArrays.Fill{
#            typeof(Gridap.Fields.affine_map), 1, Tuple{Base.OneTo{Int64}}
#        }, Gridap.Fields.AffineField{3, 3, Float64, 9}, 1, Tuple{
#            LazyArray{
#                FillArrays.Fill{
#                    Gridap.Fields.LinearCombinationMap{Colon}, 1, Tuple{Base.OneTo{Int64}}
#                }, TensorValue{3, 3, Float64, 9}, 1, Tuple{
#                    LazyArray{
#                        FillArrays.Fill{
#                            Broadcasting{
#                                Reindex{Vector{VectorValue{3, Float64}}}
#                            }, 1, Tuple{Base.OneTo{Int64}}
#                        }, Vector{VectorValue{3, Float64}}, 1, Tuple{Table{Int64, Vector{Int64}, Vector{Int32}}}
#                    }, CompressedArray{
#                        Vector{VectorValue{3, Float64}}, 1, Vector{Vector{VectorValue{3, Float64}}}, Vector{Int8}
#                    }
#                }
#            }, LazyArray{
#                FillArrays.Fill{
#                    Gridap.Fields.LinearCombinationMap{Colon}, 1, Tuple{Base.OneTo{Int64}}
#                }, VectorValue{3, Float64}, 1, Tuple{
#                    LazyArray{
#                        FillArrays.Fill{
#                            Broadcasting{
#                                Reindex{Vector{VectorValue{3, Float64}}}
#                            }, 1, Tuple{Base.OneTo{Int64}}
#                        }, Vector{VectorValue{3, Float64}}, 1, Tuple{
#                            Table{
#                                Int64, Vector{Int64}, Vector{Int32}
#                            }
#                        }
#                    }, CompressedArray{
#                        Vector{Float64}, 1, Vector{Vector{Float64}}, Vector{Int8}
#                    }
#                }
#            }
#        }
#    }, 
#    CachedVector{Int32, Vector{Int32}}
#}