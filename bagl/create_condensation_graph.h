// Copyright 2002 Indiana University.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek

#ifndef BAGL_BAGL_CREATE_CONDENSATION_GRAPH_H_
#define BAGL_BAGL_CREATE_CONDENSATION_GRAPH_H_

#include <algorithm>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"

namespace bagl {

template <concepts::AdjacencyGraph G, typename ComponentLists,
          concepts::ReadableVertexPropertyMap<G> ComponentNumberMap, concepts::EdgeMutableGraph CondGraph,
          concepts::ReadWriteEdgePropertyMap<CondGraph> EdgeMultiplicityMap>
void create_condensation_graph(const G& g, const ComponentLists& components, ComponentNumberMap component_number,
                               CondGraph& cg, EdgeMultiplicityMap edge_mult_map) {
  using CGVertex = graph_vertex_descriptor_t<CondGraph>;
  std::vector<CGVertex> to_cg_vertex(components.size());
  for (std::size_t s = 0; s < components.size(); ++s) {
    to_cg_vertex[s] = add_vertex(cg);
  }

  for (std::size_t si = 0; si < components.size(); ++si) {
    CGVertex s = to_cg_vertex[si];
    std::vector<CGVertex> adj;
    for (std::size_t i = 0; i < components[si].size(); ++i) {
      auto u = components[s][i];
      for (auto v : adjacent_vertices(u, g)) {
        CGVertex t = to_cg_vertex[component_number[v]];
        if (s != t) {
          // Avoid loops in the condensation graph
          adj.push_back(t);
        }
      }
    }
    std::sort(adj.begin(), adj.end());
    if (!adj.empty()) {
      std::size_t i = 0;
      CGVertex t = adj[i];
      auto [e, inserted] = add_edge(s, t, cg);
      put(edge_mult_map, e, 1);
      ++i;
      while (i < adj.size()) {
        if (adj[i] == t) {
          put(edge_mult_map, e, get(edge_mult_map, e) + 1);
        } else {
          t = adj[i];
          std::tie(e, inserted) = add_edge(s, t, cg);
          put(edge_mult_map, e, 1);
        }
        ++i;
      }
    }
  }
}

template <concepts::AdjacencyGraph G, typename ComponentLists,
          concepts::ReadableVertexPropertyMap<G> ComponentNumberMap, concepts::EdgeMutableGraph CondGraph>
void create_condensation_graph(const G& g, const ComponentLists& components, ComponentNumberMap component_number,
                               CondGraph& cg) {
  create_condensation_graph(g, components, component_number, cg, dummy_property_map());
}

}  // namespace bagl

#endif  // BAGL_BAGL_CREATE_CONDENSATION_GRAPH_H_
