// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PRIM_MINIMUM_SPANNING_TREE_H_
#define BAGL_BAGL_PRIM_MINIMUM_SPANNING_TREE_H_

#include <functional>

#include "bagl/graph_traits.h"
#include "bagl/dijkstra_shortest_paths.h"

namespace bagl {

template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWritePropertyMap<graph_vertex_descriptor_t<G>> PredecessorMap,
          concepts::ReadWritePropertyMap<graph_vertex_descriptor_t<G>> DistanceMap,
          concepts::ReadablePropertyMap<graph_edge_descriptor_t<G>> WeightMap,
          concepts::ReadWritePropertyMap<graph_vertex_descriptor_t<G>> IndexMap>
void prim_minimum_spanning_tree(const G& g, graph_vertex_descriptor_t<G> start, PredecessorMap predecessor,
                                DistanceMap distance, WeightMap weight, IndexMap index_map, V vis) {
  using W = property_traits_value_t<WeightMap>;
  auto compare = std::less<>{};
  auto combine = [](auto, auto b) { return b; };
  dijkstra_shortest_paths(g, start, predecessor, distance, weight, index_map, compare, combine,
                          std::numeric_limits<W>::max(), W{}, vis);
}

}  // namespace bagl

#endif // BAGL_BAGL_PRIM_MINIMUM_SPANNING_TREE_H_
