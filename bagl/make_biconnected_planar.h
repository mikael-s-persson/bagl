// Copyright 2007 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_MAKE_BICONNECTED_PLANAR_H_
#define BAGL_BAGL_MAKE_BICONNECTED_PLANAR_H_

#include <algorithm>
#include <iterator>
#include <limits>
#include <tuple>
#include <vector>

#include "bagl/biconnected_components.h"
#include "bagl/planar_detail/add_edge_visitors.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"

namespace bagl {

template <concepts::IncidenceGraph G, concepts::LvalueVertexPropertyMap<G> PlanarEmbedding,
          concepts::ReadableEdgeIndexMap<G> EdgeIndexMap, typename AddEdgeVisitor>
void make_biconnected_planar(G& g, PlanarEmbedding embedding, EdgeIndexMap em, AddEdgeVisitor& vis) {
  using vertex_t = graph_vertex_descriptor_t<G>;

  std::size_t n_edges(num_edges_or_zero(g));
  std::vector<vertex_t> articulation_points;
  auto component_map = vector_property_map(n_edges, em, std::size_t{0});

  biconnected_components(g, component_map, std::back_inserter(articulation_points));

  for (auto v : articulation_points) {
    std::size_t previous_component = std::numeric_limits<std::size_t>::max();
    vertex_t previous_vertex = graph_traits<G>::null_vertex();

    for (auto e : embedding[v]) {
      vertex_t e_source(source(e, g));
      vertex_t e_target(target(e, g));

      // Skip self-loops and parallel edges
      if (e_source == e_target || previous_vertex == e_target) {
        continue;
      }

      vertex_t current_vertex = e_source == v ? e_target : e_source;
      std::size_t current_component = component_map[e];
      if (previous_vertex != graph_traits<G>::null_vertex() && current_component != previous_component) {
        vis.visit_vertex_pair(current_vertex, previous_vertex, g);
      }
      previous_vertex = current_vertex;
      previous_component = current_component;
    }
  }
}

template <concepts::IncidenceGraph G, concepts::LvalueVertexPropertyMap<G> PlanarEmbedding,
          concepts::ReadableEdgeIndexMap<G> EdgeIndexMap>
void make_biconnected_planar(G& g, PlanarEmbedding embedding, EdgeIndexMap em) {
  default_add_edge_visitor vis;
  make_biconnected_planar(g, embedding, em, vis);
}

template <concepts::IncidenceGraph G, concepts::LvalueVertexPropertyMap<G> PlanarEmbedding>
void make_biconnected_planar(G& g, PlanarEmbedding embedding) {
  make_biconnected_planar(g, embedding, get(edge_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_MAKE_BICONNECTED_PLANAR_H_
