// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PRIM_MINIMUM_SPANNING_TREE_H_
#define BAGL_BAGL_PRIM_MINIMUM_SPANNING_TREE_H_

#include <functional>

#include "bagl/dijkstra_shortest_paths.h"
#include "bagl/graph_traits.h"
#include "bagl/vector_property_map.h"

namespace bagl {

template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> IndexMap>
void prim_minimum_spanning_tree(const G& g, graph_vertex_descriptor_t<G> start, PredecessorMap predecessor,
                                DistanceMap distance, WeightMap weight, IndexMap index_map, V vis) {
  using W = property_traits_value_t<WeightMap>;
  auto compare = std::less<>{};
  auto combine = [](auto, auto b) { return b; };
  dijkstra_shortest_paths(g, start, predecessor, distance, weight, index_map, compare, combine,
                          std::numeric_limits<W>::max(), W{}, vis);
}

template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, concepts::ReadWriteVertexPropertyMap<G> IndexMap>
void prim_minimum_spanning_tree(const G& g, graph_vertex_descriptor_t<G> start, PredecessorMap predecessor,
                                WeightMap weight, IndexMap index_map) {
  using W = property_traits_value_t<WeightMap>;
  auto compare = std::less<>{};
  auto combine = [](auto, auto b) { return b; };
  dijkstra_shortest_paths(g, start, predecessor, vector_property_map(num_vertices(g), index_map, W{}).ref(), weight,
                          index_map, compare, combine, std::numeric_limits<W>::max(), W{}, default_dijkstra_visitor());
}

}  // namespace bagl

#endif // BAGL_BAGL_PRIM_MINIMUM_SPANNING_TREE_H_
